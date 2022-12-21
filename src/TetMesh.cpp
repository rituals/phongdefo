#include "TetMesh.h"

#include <nanoflann.hpp>
#include <Eigen/Dense>
#include <array>
#include <set>
#include <vector>
#include <functional>

using namespace TetEmbeddedDeform;

bool
Embedding::isInside() const
{
    return (
        coords[0] >= 0.0f && coords[0] <= 1.0f &&
        coords[1] >= 0.0f && coords[1] <= 1.0f &&
        coords[2] >= 0.0f && coords[2] <= 1.0f &&
        coords[3] >= 0.0f && coords[3] <= 1.0f);
}

TetMesh::TetMesh() :
    _points(),
    _cellVerts(),
    _vertCells(),
    _centroids(),
    _edgeBases(),
    _edgeBasesInv()
{
}

TetMesh::~TetMesh()
{
}

void
TetMesh::setPoints(std::vector<Vector3> &&points)
{
    _points = std::move(points);
}

void
TetMesh::setVerts(std::vector<std::array<size_t, 4>> &&verts)
{
    _cellVerts = std::move(verts);
}

void
TetMesh::build()
{
    computeVertCells();
    computeCentroids();
    computeEdgeBases(true);
}

size_t
TetMesh::getNumPoints() const
{
    return _points.size();
}

size_t
TetMesh::getNumCells() const
{
    return _cellVerts.size(); 
}

const std::vector<Vector3> &
TetMesh::getPoints() const
{
    return _points;
}

const Vector3 &
TetMesh::getPoint(size_t index) const
{
    return _points[index];
}

const std::array<size_t, 4> &
TetMesh::getCellVerts(size_t cell) const
{
    return _cellVerts[cell];
}

const std::vector<size_t> &
TetMesh::getVertCells(size_t vert) const
{
    return _vertCells[vert];
}

Vector3
TetMesh::getCentroid(size_t cell) const
{
    return _centroids[cell];
}

Matrix3
TetMesh::getEdgeBasis(size_t cell) const
{
    return _edgeBases[cell];
}

Matrix3
TetMesh::getEdgeBasisInv(size_t cell) const
{
    return _edgeBasesInv[cell];
}

std::vector<Embedding>
TetMesh::computeEmbedding(const std::vector<Vector3> &points) const
{
    Matrix3X pointsMat;
    pointsMat.resize(3, _points.size());
    for (size_t i = 0; i < _points.size(); ++i) {
        pointsMat.col(i) = _points[i];
    }

    // Create kdtree for faster searching through the tet mesh.
    using kdtree = nanoflann::KDTreeEigenMatrixAdaptor<
        Matrix3X, 3, nanoflann::metric_L2, false>;

    kdtree index(3, std::cref(pointsMat), 5);

    std::vector<Embedding> embeddings(points.size());

    for (size_t i = 0; i < points.size(); ++i) {

        // Query the nearest cell vertices. 
        const size_t nresults = 4;
        std::vector<size_t> indices(nresults);
        std::vector<float> dists(nresults);
        nanoflann::KNNResultSet<float> results(nresults);
        results.init(&indices[0], &dists[0]);
        index.index_->findNeighbors(results, points[i].data());

        // Map those points to coincident cells. 
        std::vector<size_t> cells;
        std::set<size_t> seen;
        for (size_t m = 0; m < nresults; ++m) {
            //for (size_t cell : _vertCells[matches[m].first]) {
            for (size_t cell : _vertCells[indices[m]]) {
                // Matches are ordered by proximity so preserve that order while
                // removing dupes.
                if (seen.find(cell) == seen.end()) {
                    seen.insert(cell);
                    cells.push_back(cell);
                }
            }
        }

        // Compute the barycentric coords with respect to the gathered cells.
        // Stop when the point is found inside a cell.
        for (size_t cell : cells) {
            embeddings[i].cell = cell;
            computeBarycentricCoords(points[i], &embeddings[i]);
            if (embeddings[i].isInside()) {
                break;
            }
        }

        // Check if no enclosing cell was found. This can happen for points
        // exterior to the tet mesh. Anchor these points to the first matched
        // cell. This seems to work ok in practice but something to keep an eye
        // on.
        if (!embeddings[i].isInside()) {
            embeddings[i].cell = cells[0];
            computeBarycentricCoords(points[i], &embeddings[i]);
        }
    }

    return embeddings;
}

void
TetMesh::computeVertCells()
{
    _vertCells.assign(_points.size(), std::vector<size_t>());

    for (size_t cell = 0; cell < _cellVerts.size(); ++cell) {
        for (auto vert : _cellVerts[cell]) {
            _vertCells[vert].push_back(cell);
        }
    }
}

void
TetMesh::computeCentroids()
{
    const size_t ncells = _cellVerts.size();
    _centroids.resize(ncells);

    for (size_t cell = 0; cell < ncells; ++cell) {
        _centroids[cell] = (
            _points[_cellVerts[cell][0]] +
            _points[_cellVerts[cell][1]] +
            _points[_cellVerts[cell][2]] +
            _points[_cellVerts[cell][3]]) / 4.0;
    }
}

void
TetMesh::computeEdgeBases(bool computeInv)
{
    const size_t ncells = _cellVerts.size();
    _edgeBases.resize(ncells);
    _edgeBasesInv.resize(ncells);

    for (size_t cell = 0; cell < _cellVerts.size(); ++cell) {
        Matrix3 m;
        m << _points[_cellVerts[cell][1]] - _points[_cellVerts[cell][0]],
             _points[_cellVerts[cell][2]] - _points[_cellVerts[cell][0]],
             _points[_cellVerts[cell][3]] - _points[_cellVerts[cell][0]];

        _edgeBases[cell] = m;
        if (computeInv) {
            _edgeBasesInv[cell] = m.inverse();
        } else {
            _edgeBasesInv[cell] = Matrix3::Identity();
        }
    }
}

void
TetMesh::computeBarycentricCoords(const Vector3 &point, Embedding *e) const
{
    // Refer to eq. (4), (5).
    const auto &cellVerts = _cellVerts[e->cell];
    const Vector3 x0 = _points[cellVerts[0]];
    const Vector3 coords = _edgeBasesInv[e->cell]*(point - x0);
    e->coords[0] = 1.0f - coords[0] - coords[1] - coords[2];
    e->coords[1] = coords[0];
    e->coords[2] = coords[1];
    e->coords[3] = coords[2];
}
