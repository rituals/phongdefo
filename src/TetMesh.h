#ifndef __TETMESH_H__
#define __TETMESH_H__

#include <Eigen/Dense>

namespace TetEmbeddedDeform
{

using Vector3 = Eigen::Vector3f;
using Matrix3 = Eigen::Matrix3f;
using Matrix3X = Eigen::Matrix3Xf; 

// The embedding of a point inside a tetrahedral mesh.
struct Embedding
{
    size_t cell;        // index of the enclosing cell in the tet mesh.
    float coords[4];    // barycentric coords.

    bool isInside() const;
};

// A tetrahedral mesh, defined simply via points and 4-tuple cell vertices.
class TetMesh
{
public:
    TetMesh();
    ~TetMesh();

    void setPoints(std::vector<Vector3> &&points);

    void setVerts(std::vector<std::array<size_t, 4>> &&verts);

    // Call this before querying anything.
    void build();

    size_t getNumPoints() const;

    size_t getNumCells() const;

    const std::vector<Vector3> &getPoints() const;

    const Vector3 &getPoint(size_t index) const;

    // Returns the 4 corner vertices of the cell.
    const std::array<size_t, 4> &getCellVerts(size_t cell) const;

    // Returns the cells adjacent to the vertex.
    const std::vector<size_t> &getVertCells(size_t vert) const;

    // Returns the centroid of the cell.
    Vector3 getCentroid(size_t cell) const;

    // Returns the 3x3 edge basis matrix of the cell.
    Matrix3 getEdgeBasis(size_t cell) const;

    // Returns the inverted 3x3 edge basis matrix of the cell.
    Matrix3 getEdgeBasisInv(size_t cell) const;

    // Computes and returns the embedding for the points.
    std::vector<Embedding> computeEmbedding(
        const std::vector<Vector3> &points) const;

    void computeVertCells();

    void computeCentroids();

    void computeEdgeBases(bool computeInv);

    void computeBarycentricCoords(const Vector3 &point, Embedding *e) const;

private:
    // Point positions.
    std::vector<Vector3> _points;

    // Mapping from cell to its 4 corner vertices.
    std::vector<std::array<size_t, 4>> _cellVerts;

    // Mapping from vertex to coincident cells.
    std::vector<std::vector<size_t>> _vertCells;

    // Cell centroids.
    std::vector<Vector3> _centroids;

    // Cell edge bases. 
    std::vector<Matrix3> _edgeBases;

    // Inverse cell edge bases.
    std::vector<Matrix3> _edgeBasesInv;
};

} // end namespace TetEmbeddedDeform

#endif
