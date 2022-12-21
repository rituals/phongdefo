#include "DeformationTracker.h"

#include <UT/UT_Assert.h>
#include <array>
#include <vector>
#include <iostream>

using namespace TetEmbeddedDeform;

static std::vector<Vector3> 
extractPoints(const GU_Detail *detail);

static std::vector<std::array<size_t, 4>>
extractCellVerts(const GU_Detail *detail); 

DeformationTracker::DeformationTracker() :
    _restMeshDetail(),
    _deformedMeshDetail(),
    _restMesh(),
    _deformedMesh(),
    _points(),
    _embedding(),
    _vertDGs(),
    _cellDGs(),
    _vertDGweights(),
    _restPointsId(-1),
    _restTopoId(-1),
    _defPointsId(-1),
    _defTopoId(-1)
{
}

DeformationTracker::~DeformationTracker()
{
}

void
DeformationTracker::update(
    const GU_Detail *detail,
    const GU_Detail *restMeshDetail,
    const GU_Detail *deformedMeshDetail)
{
    _restMeshDetail = restMeshDetail;
    _deformedMeshDetail = deformedMeshDetail;

    auto restPointsId = restMeshDetail->getP()->getDataId();
    auto restTopoId = restMeshDetail->getTopology().getDataId();
    auto defPointsId = deformedMeshDetail->getP()->getDataId();
    auto defTopoId = deformedMeshDetail->getTopology().getDataId();

    bool rebuildRestMesh =
        restPointsId != _restPointsId || restTopoId != _restTopoId;

    bool rebuildDeformedMesh =
        rebuildRestMesh || defTopoId != _defTopoId;

    if (rebuildDeformedMesh) {
        buildDeformedMesh();
    }

    if (rebuildRestMesh) {
        buildRestMesh();
        _embedding = _restMesh.computeEmbedding(extractPoints(detail));
        computeVertDGweights();
    }

    computeDeformation();

    _restPointsId = restPointsId;
    _restTopoId = restTopoId;
    _defPointsId = defPointsId;
    _defTopoId = defTopoId;
}

void
DeformationTracker::computeDeformation()
{
    _deformedMesh.setPoints(extractPoints(_deformedMeshDetail));
    _deformedMesh.computeCentroids();
    _deformedMesh.computeEdgeBases(false);

    computeCellDGs(); 
    computeVertDGs();
}

UT_Vector3
DeformationTracker::applyLinearDeformation(
    const UT_Vector3 &point, size_t index) const
{
    // eq (7): f_c(X) = c + F_c(X - C)
    //
    // f_c: linear interpolation evaluated at cell center.
    // X: point to deform.
    // c: deformed cell centroid.
    // C: undeformed cell centroid.
    // F_c: cell deformation gradient.
    Vector3 X(point[0], point[1], point[2]);
    const Embedding &e = _embedding[index];
    const Vector3 &c = _deformedMesh.getCentroid(e.cell);
    const Vector3 &C = _restMesh.getCentroid(e.cell);
    const Matrix3 &F_C = _cellDGs[e.cell];
    const Vector3 f_C = c + F_C*(X - C);

    return UT_Vector3(f_C[0], f_C[1], f_C[2]);
}

UT_Vector3
DeformationTracker::applyPhongDeformation(
    const UT_Vector3 &point, size_t index) const
{
    UT_Vector3 linear = applyLinearDeformation(point, index);

    // eq (10): f_V(X) = sum(B_i*(x_i + F_i*(X - X_i))
    //
    // f_V: interpolated deformers at cell vertices.
    // X_i: rest position of vertex.
    // x_i: deformed position of vertex.
    // F_i: deformation gradient at vertex.
    Vector3 X(point[0], point[1], point[2]);
    Vector3 f_V = Vector3::Zero();
    const Embedding &e = _embedding[index];
    const auto cellVerts = _restMesh.getCellVerts(e.cell);
    for (size_t i = 0; i < 4; ++i) {
        const size_t vert = cellVerts[i];
        const Vector3 &X_i = _restMesh.getPoint(vert);
        const Vector3 &x_i = _deformedMesh.getPoint(vert);
        const Matrix3 &F_i = _vertDGs[vert];
        f_V += e.coords[i]*(x_i + F_i*(X - X_i));
    }

    // eq (12): f_phong(X) = (f_C(X) + f_V(X))/2
    return (linear + UT_Vector3(f_V[0], f_V[1], f_V[2])) / 2.0;
}

void
DeformationTracker::computeCellDGs()
{
    const size_t ncells = _restMesh.getNumCells(); 
    _cellDGs.resize(ncells);

    for (size_t cell = 0; cell < ncells; ++cell) {
        _cellDGs[cell] =
            _deformedMesh.getEdgeBasis(cell)*_restMesh.getEdgeBasisInv(cell);
    }
}

void
DeformationTracker::computeVertDGs()
{
    const size_t nverts = _restMesh.getNumPoints();
    _vertDGs.resize(nverts);

    for (size_t vert = 0; vert < nverts; ++vert) {
        _vertDGs[vert] = Matrix3::Zero();
        const auto &cells = _restMesh.getVertCells(vert);
        for (size_t i = 0; i < cells.size(); ++i) {
            _vertDGs[vert] += _vertDGweights[vert][i] * _cellDGs[cells[i]];
        }
    }
}

void
DeformationTracker::computeVertDGweights()
{
    const size_t nverts = _restMesh.getNumPoints();
    _vertDGweights.assign(nverts, std::vector<float>());

    for (size_t vert = 0; vert < nverts; ++vert) {
        const Vector3 &point = _restMesh.getPoint(vert);        
        const std::vector<size_t> &cells = _restMesh.getVertCells(vert);

        std::vector<Vector3> point_to_centroids(cells.size());
        Matrix3 A = Matrix3::Zero();
        Vector3 b = Vector3::Zero();

        for (size_t k = 0; k < cells.size(); ++k) {
            size_t cell = cells[k]; 
            point_to_centroids[k] = _restMesh.getCentroid(cell) - point;
            Vector3 p2c_norm = point_to_centroids[k].normalized();
            A += p2c_norm * p2c_norm.transpose();
            b -= p2c_norm;
        }

        A += Matrix3::Identity();
        Vector3 lambda = A.inverse()*b;

        float sum = 0.0f;
        for (size_t k = 0; k < cells.size(); ++k) {
            float len = point_to_centroids[k].norm();
            float w = 1.0f + lambda.dot(point_to_centroids[k]/len);
            w /= len;
            sum += w;
            _vertDGweights[vert].push_back(w);
        }

        for (auto &w : _vertDGweights[vert]) {
            w /= sum;
        }
    }
}

void
DeformationTracker::buildRestMesh()
{
    _restMesh.setPoints(extractPoints(_restMeshDetail));
    _restMesh.setVerts(extractCellVerts(_restMeshDetail));
    _restMesh.build();
}

void
DeformationTracker::buildDeformedMesh()
{
    _deformedMesh.setPoints(extractPoints(_deformedMeshDetail));
    _deformedMesh.setVerts(extractCellVerts(_deformedMeshDetail));
    _deformedMesh.build();
}

std::vector<Vector3>
extractPoints(const GU_Detail *detail)
{
    const GA_Size numPoints = detail->getNumPoints();
    std::vector<Vector3> points;
    points.reserve(numPoints);
    GA_ROHandleV3 pHandle(detail->getP());

    for (GA_Size i = 0; i < numPoints; ++i) {
        const auto p = pHandle(detail->pointOffset(i));
        points.push_back(Vector3(p[0], p[1], p[2]));
    }
    
    return points;
}

std::vector<std::array<size_t, 4>>
extractCellVerts(const GU_Detail *detail)
{
    std::vector<std::array<size_t, 4>> verts;
    verts.reserve(detail->getNumPrimitives());
    const GEO_Primitive *prim;

    GA_FOR_ALL_PRIMITIVES(detail, prim) {
        std::array<size_t, 4> primVerts;

        for (GA_Size i = 0; i < prim->getVertexCount(); ++i) {
            const GA_Offset vOffset = prim->getVertexOffset(i);
            const GA_Offset pOffset = detail->vertexPoint(vOffset);
            const size_t pIndex = detail->pointIndex(pOffset);
            primVerts[i] = pIndex;
        }

        verts.push_back(primVerts);
    }
    
    return verts;
}
