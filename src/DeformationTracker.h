#ifndef __DEFORMATIONTRACKER_H__
#define __DEFORMATIONTRACKER_H__

#include "TetMesh.h"
#include <GA/GA_Types.h>
#include <GU/GU_Detail.h>
#include <UT/UT_Vector3.h>
#include <vector>

namespace TetEmbeddedDeform
{

// Class for tracking tetrahedral mesh deformation and applying it to a set
// of interior points. The deformation is defined by the delta of two 
// tetrahedral meshes, rest and deformed, whose topology must match.
class DeformationTracker
{
public:
    DeformationTracker();

    ~DeformationTracker();

    void update(
        const GU_Detail *detail,
        const GU_Detail *restMeshDetail,
        const GU_Detail *deformedMeshDetail);

    UT_Vector3 applyLinearDeformation(
        const UT_Vector3 &point, size_t index) const;

    UT_Vector3 applyPhongDeformation(
        const UT_Vector3 &point, size_t index) const;

private:
    // Computes the cell-based deformation gradients.
    void computeCellDGs();

    // Computes the vertex-based deformation gradients.
    void computeVertDGs();

    // Computes the weights used for interpolating vertex deformation
    // gradients from incident cells.
    void computeVertDGweights();

    void buildRestMesh();

    void buildDeformedMesh();

    void computeDeformation();

private:
    const GU_Detail *_restMeshDetail;
    const GU_Detail *_deformedMeshDetail;

    TetMesh _restMesh;
    TetMesh _deformedMesh;

    // Interior points.
    std::vector<Vector3> _points;

    // Embedding of interior points.
    std::vector<Embedding> _embedding;

    // Vertex-based deformation gradients.
    std::vector<Matrix3> _vertDGs;

    // Cell-based deformation gradients.
    std::vector<Matrix3> _cellDGs;

    // Weights used for interpolating vertex DGs from cell DGs.
    std::vector<std::vector<float>> _vertDGweights;

    // For tracking changes on the rest, deformed tet meshes.
    GA_DataId _restPointsId;
    GA_DataId _restTopoId;
    GA_DataId _defPointsId;
    GA_DataId _defTopoId;
};

} // end namespace TetEmbeddedDeform

#endif
