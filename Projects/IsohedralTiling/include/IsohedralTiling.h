#ifndef ISOHEDRAL_TILING_H
#define ISOHEDRAL_TILING_H

#include <utility>
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <tbb/tbb.h>
#include <unordered_map>
#include <unordered_set>
#include <complex>
#include <iomanip>

#include <random>
#include <cmath>
#include <fstream>
#include <gmsh.h>

#include "../tactile/tiling.hpp"
#include "../clipper/clipper.hpp"

#include "VecMatDef.h"

class IsohedralTiling
{
public:
    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
    using MatrixXT = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXi = Vector<int, Eigen::Dynamic>;
    using MatrixXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;

    using TV = Vector<T, 3>;
    using TV2 = Vector<T, 2>;
    using TM2 = Matrix<T, 2, 2>;
    using TV3 = Vector<T, 3>;
    using IV = Vector<int, 3>;
    using IV2 = Vector<int, 2>;
    using TM = Matrix<T, 3, 3>;
    using StiffnessMatrix = Eigen::SparseMatrix<T>;
    using Entry = Eigen::Triplet<T>;


    using dvec2 = glm::dvec2;
    using dmat3 = glm::dmat3;

public:
    VectorXT deformed;
    VectorXi faces;

    std::string data_folder = "./";

    bool quadratic = false;
public:
    IsohedralTiling() {}
    ~IsohedralTiling() {}

    bool fetchUnitCellFromOneFamily(int IH, int n_unit,
        std::vector<std::vector<TV2>>& eigen_polygons,
        std::vector<TV2>& eigen_base, 
        const std::vector<T>& params,
        const Vector<T, 4>& eij, const std::string& filename,
        T unit);
    void generatePeriodicMesh(std::vector<std::vector<TV2>>& polygons, 
        std::vector<TV2>& pbc_corners, bool save_to_file, std::string prefix);
    void generateOnePerodicUnit();

private:
    T evalDistance(const TV& p1, const TV& p2, const TV& q, T t)
    {
        return std::sqrt(std::pow((0.1e1 - t) * p1[0] + t * p2[0] - q[0], 0.2e1) + std::pow((0.1e1 - t) * p1[1] + t * p2[1] - q[1], 0.2e1));
    };

    T closestTToLine(const TV& p1, const TV& p2, const TV& q){
        return (p1[0] * p1[0] + (-p2[0] - q[0]) * p1[0] + p1[1] * p1[1] + (-p2[1] - q[1]) * p1[1] + p2[0] * q[0] + p2[1] * q[1]) / (p1[0] * p1[0] - 2 * p1[0] * p2[0] + p1[1] * p1[1] - 2 * p1[1] * p2[1] + p2[0] * p2[0] + p2[1] * p2[1]);
    };
    // operations on the tiling unit
    glm::dmat3 centrePSRect(T xmin, T ymin, T xmax, T ymax);
    std::vector<glm::dvec2> outShapeVec(const std::vector<glm::dvec2>& vec, const glm::dmat3& M);
    void getTilingShape(std::vector<dvec2>& shape, const csk::IsohedralTiling& tiling,
        const std::vector<std::vector<dvec2>>& edges);
    void getTilingEdges(const csk::IsohedralTiling& tiling,
        const Vector<T, 4>& eij,
        std::vector<std::vector<dvec2>>& edges);
    void getTranslationUnitPolygon(std::vector<std::vector<dvec2>>& polygons_v,
        const std::vector<dvec2>& shape,
        const csk::IsohedralTiling& tiling, Vector<T, 4>& transf,
        int width, int depth, TV2& xy_shift);
    void shapeToPolygon(ClipperLib::Paths& final_shape, 
        std::vector<std::vector<TV2>>& polygons, T mult);
    void periodicToBase(const Vector<T, 8>& periodic, std::vector<TV2>& eigen_base);
    void saveClip(const ClipperLib::Paths& final_shape, 
        const Vector<T, 8>& periodic, T mult,
        const std::string& filename, bool add_box = true);
};

#endif