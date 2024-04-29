#include "../include/IsohedralTiling.h"
#include "../include/Util.h"

// ============================ PUBLIC ============================
void IsohedralTiling::generateOnePerodicUnit()
{
    std::vector<std::vector<TV2>> polygons;
    std::vector<TV2> pbc_corners; 
    int tiling_idx = 0;

    T unit = 5.0;
    if (tiling_idx == 0 || tiling_idx == 19)
        unit = 5.0;
    else if (tiling_idx == 26)
        unit = 6.0;
    else if (tiling_idx == 27 )
        unit = 10.0;
    else if (tiling_idx == 60)
        unit = 6.0;
    else if (tiling_idx == 46)
        unit = 3.0;
    else if (tiling_idx == 20)
        unit = 7.0;
    csk::IsohedralTiling a_tiling( csk::tiling_types[ tiling_idx ] );
    int num_params = a_tiling.numParameters();
    T new_params[ num_params ];
    a_tiling.getParameters( new_params );
    std::vector<T> params(num_params);
    for (int j = 0; j < num_params;j ++)
        params[j] = new_params[j];
    
    Vector<T, 4> cubic_weights;
    cubic_weights << 0.25, 0., 0.75, 0.;
    bool non_ixn_unit = fetchUnitCellFromOneFamily(tiling_idx, 2, polygons, pbc_corners, params, 
        cubic_weights, data_folder + "a_structure.txt", unit);
    if (!non_ixn_unit)
    {
        std::cout << "self intersecting unit" << std::endl;
        return;
    }

    generatePeriodicMesh(polygons, pbc_corners, true, data_folder + "a_structure");

    Eigen::MatrixXd V; Eigen::MatrixXi F, V_quad;

    if (quadratic)
    {
        loadQuadraticTriangleMeshFromVTKFile(data_folder + "a_structure.vtk", V, F, V_quad);
        F.resize(V_quad.rows(), 3);
        F.col(0) = V_quad.col(0); F.col(1) = V_quad.col(1); F.col(2) = V_quad.col(2);
        TV3 e0(V.row(F(0, 1)) - V.row(F(0, 0)));
        TV3 e1(V.row(F(0, 2)) - V.row(F(0, 0)));
        if (e1.cross(e0).dot(TV3(0, 0, 1)) > 0)
        {
            F.col(0) = V_quad.col(0); F.col(1) = V_quad.col(2); F.col(2) = V_quad.col(1);
            Eigen::MatrixXi V_quad_backup = V_quad;
            V_quad.col(1) = V_quad_backup.col(2); V_quad.col(2) = V_quad_backup.col(1);
            V_quad.col(5) = V_quad_backup.col(4); V_quad.col(4) = V_quad_backup.col(5);
        }
    }
    else
    {
        loadMeshFromVTKFile(data_folder + "a_structure.vtk", V, F);
        iglMatrixFatten<T, 3>(V, deformed);
        iglMatrixFatten<int, 3>(F, faces);
    }
}

void IsohedralTiling::generatePeriodicMesh(std::vector<std::vector<TV2>>& polygons, 
    std::vector<TV2>& pbc_corners, bool save_to_file, std::string prefix)
{
      // Before using any functions in the C++ API, Gmsh must be initialized:

    TV2 p1 = pbc_corners[0];
    TV2 p2 = pbc_corners[1];
    TV2 p3 = pbc_corners[2];
    TV2 p4 = pbc_corners[3];

    T eps = 1e-6;
    gmsh::initialize();

    gmsh::model::add("tiling");
    // gmsh::logger::start();
    // gmsh::logger::stop();

    gmsh::option::setNumber("Geometry.ToleranceBoolean", eps);
    gmsh::option::setNumber("Geometry.Tolerance", eps);
    if (quadratic)
        gmsh::option::setNumber("Mesh.ElementOrder", 2);
    else
        gmsh::option::setNumber("Mesh.ElementOrder", 1);

    gmsh::option::setNumber("General.Verbosity", 0);

    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    T th = eps;
    //Points
    int acc = 1;

    // clamping box 1 2 3 4
    for (int i = 0; i < pbc_corners.size(); ++i)
		gmsh::model::occ::addPoint(pbc_corners[i][0], pbc_corners[i][1], 0, 1, acc++);

    for(int i=0; i<polygons.size(); ++i)
    {
        for(int j=0; j<polygons[i].size(); ++j)
        {
            gmsh::model::occ::addPoint(polygons[i][j][0], polygons[i][j][1], 0, 1, acc++);
        }
    }
    
    //Lines
    acc = 1;
    int starting_vtx = 1;

    // add clipping box
    gmsh::model::occ::addLine(1, 2, acc++); 
    gmsh::model::occ::addLine(2, 3, acc++); 
    gmsh::model::occ::addLine(3, 4, acc++); 
    gmsh::model::occ::addLine(4, 1, acc++);

    starting_vtx = 5;
    std::vector<T> poly_lines;
    for (int i = 0; i < polygons.size(); i++)
    {
        for(int j=1; j<polygons[i].size(); ++j)
        {
            gmsh::model::occ::addLine(starting_vtx++, starting_vtx, acc++);
            poly_lines.push_back(acc);
        }
        gmsh::model::occ::addLine(starting_vtx, starting_vtx-polygons[i].size()+1, acc++);
        poly_lines.push_back(acc);
        ++starting_vtx;
    }

    gmsh::model::mesh::field::add("Distance", 1);
    // gmsh::model::mesh::field::setNumbers(1, "CurvesList", poly_lines);

    gmsh::model::mesh::field::add("Threshold", 2);
    gmsh::model::mesh::field::setNumber(2, "InField", 1);
    gmsh::model::mesh::field::setNumber(2, "SizeMin", 0.2);
    gmsh::model::mesh::field::setNumber(2, "SizeMax", 1.0);
    gmsh::model::mesh::field::setNumber(2, "DistMin", 0.005);

    // gmsh::model::mesh::field::setNumber(2, "SizeMin", 1.0);
    // gmsh::model::mesh::field::setNumber(2, "SizeMax", 2.0);
    // gmsh::model::mesh::field::setNumber(2, "DistMin", 0.005);

    gmsh::model::mesh::field::setAsBackgroundMesh(2);
    
    acc = 1;
    int acc_loop = 1;

    gmsh::model::occ::addCurveLoop({1, 2, 3, 4}, acc++);
    acc_loop = 5;
    for (int i = 0; i < polygons.size(); i++)
    {
        std::vector<int> polygon_loop;
        for(int j=1; j < polygons[i].size()+1; j++)
            polygon_loop.push_back(acc_loop++);
        gmsh::model::occ::addCurveLoop(polygon_loop, acc++);
    }

    for (int i = 0; i < polygons.size()+1; i++)
    {
        gmsh::model::occ::addPlaneSurface({i+1});
    }

    std::vector<std::pair<int ,int>> poly_idx;
    for(int i=0; i<polygons.size(); ++i)
        poly_idx.push_back(std::make_pair(2, i+2));
    
    std::vector<std::pair<int, int>> ov;
    std::vector<std::vector<std::pair<int, int> > > ovv;
    gmsh::model::occ::cut({{2, 1}}, poly_idx, ov, ovv);
    gmsh::model::occ::synchronize();

    int zero_idx;
    for(int i=0; i < pbc_corners.size(); i++)
    {
        if(pbc_corners[i].norm()<1e-6)
        {
            zero_idx = i;
            break;
        }
    }

    TV2 t1 = pbc_corners[(zero_idx+1)%pbc_corners.size()];
    TV2 t2 = pbc_corners[(zero_idx+3)%pbc_corners.size()];

    std::vector<T> translation_hor({1, 0, 0, t1[0], 0, 1, 0, t1[1], 0, 0, 1, 0, 0, 0, 0, 1});
	std::vector<T> translation_ver({1, 0, 0, t2[0], 0, 1, 0, t2[1], 0, 0, 1, 0, 0, 0, 0, 1});

    std::vector<std::pair<int, int>> sleft;
    gmsh::model::getEntitiesInBoundingBox(std::min(0.0,t2[0])-eps, std::min(0.0,t2[1])-eps, -eps, std::max(0.0,t2[0])+eps, std::max(0.0,t2[1])+eps, eps, sleft, 1);
    // std::ofstream pbc_output(data_folder + "pbc_data.txt");
    for(auto i : sleft) {
        T xmin, ymin, zmin, xmax, ymax, zmax;
        gmsh::model::getBoundingBox(i.first, i.second, xmin, ymin, zmin, xmax, ymax, zmax);
        std::vector<std::pair<int, int> > sright;
        gmsh::model::getEntitiesInBoundingBox(xmin-eps+t1[0], ymin-eps+t1[1], zmin - eps, xmax+eps+t1[0], ymax+eps+t1[1], zmax + eps, sright, 1);

        for(auto j : sright) {
            T xmin2, ymin2, zmin2, xmax2, ymax2, zmax2;
            gmsh::model::getBoundingBox(j.first, j.second, xmin2, ymin2, zmin2, xmax2, ymax2, zmax2);
            xmin2 -= t1[0];
            ymin2 -= t1[1];
            xmax2 -= t1[0];
            ymax2 -= t1[1];
            if(std::abs(xmin2 - xmin) < eps && std::abs(xmax2 - xmax) < eps &&
                std::abs(ymin2 - ymin) < eps && std::abs(ymax2 - ymax) < eps &&
                std::abs(zmin2 - zmin) < eps && std::abs(zmax2 - zmax) < eps) 
            {
                gmsh::model::mesh::setPeriodic(1, {j.second}, {i.second}, translation_hor);
                // pbc_output << "X " << j.second << " " << i.second << std::endl;
            }
        }
    }

    std::vector<std::pair<int, int>> sbottom;
    gmsh::model::getEntitiesInBoundingBox(std::min(0.0,t1[0])-eps, std::min(0.0,t1[1])-eps, -eps, std::max(0.0,t1[0])+eps, std::max(0.0,t1[1])+eps, eps, sbottom, 1);

    for(auto i : sbottom) {
        T xmin, ymin, zmin, xmax, ymax, zmax;
        gmsh::model::getBoundingBox(i.first, i.second, xmin, ymin, zmin, xmax, ymax, zmax);
        std::vector<std::pair<int, int> > stop;
        gmsh::model::getEntitiesInBoundingBox(xmin-eps+t2[0], ymin-eps+t2[1], zmin - eps, xmax+eps+t2[0], ymax+eps+t2[1], zmax + eps, stop, 1);

        for(auto j : stop) {
            T xmin2, ymin2, zmin2, xmax2, ymax2, zmax2;
            gmsh::model::getBoundingBox(j.first, j.second, xmin2, ymin2, zmin2, xmax2, ymax2, zmax2);
            xmin2 -= t2[0];
            ymin2 -= t2[1];
            xmax2 -= t2[0];
            ymax2 -= t2[1];
            if(std::abs(xmin2 - xmin) < eps && std::abs(xmax2 - xmax) < eps &&
                std::abs(ymin2 - ymin) < eps && std::abs(ymax2 - ymax) < eps &&
                std::abs(zmin2 - zmin) < eps && std::abs(zmax2 - zmax) < eps) 
            {
                gmsh::model::mesh::setPeriodic(1, {j.second}, {i.second}, translation_ver);
                // pbc_output << "Y " << j.second << " " << i.second << std::endl;
            }
        }
    }
    
    gmsh::model::occ::synchronize();

    gmsh::model::mesh::generate(2);

    
    gmsh::write(prefix + ".vtk");
    std::ofstream translation(prefix + "_translation.txt");
    translation << std::setprecision(20) << t1.transpose() << std::endl;
    translation << std::setprecision(20) << t2.transpose() << std::endl;
    translation.close();
    gmsh::finalize();
    
}

bool IsohedralTiling::fetchUnitCellFromOneFamily(int IH, int n_unit,
    std::vector<std::vector<TV2>>& eigen_polygons,
    std::vector<TV2>& eigen_base, 
    const std::vector<T>& params,
    const Vector<T, 4>& eij, const std::string& filename,
    T unit)
{
    using namespace csk;
    using namespace std;
    using namespace glm;

    std::ofstream out(filename);
    csk::IsohedralTiling a_tiling( csk::tiling_types[ IH ] );
    out << "IH " << IH << std::endl;
    size_t num_params = a_tiling.numParameters();
    T new_params[ num_params ];
    a_tiling.getParameters( new_params );
    // Change a parameter
    for( size_t idx = 0; idx < a_tiling.numParameters(); ++idx ) 
    {
        new_params[idx] = params[idx];
        out << params[idx] << " ";
    }
    out << std::endl;
    a_tiling.setParameters( new_params );

    out << "numEdgeShapes " << int(a_tiling.numEdgeShapes()) << " ";
    std::vector<std::vector<dvec2>> edges(a_tiling.numEdgeShapes());
    getTilingEdges(a_tiling, eij, edges);
    for (auto ej : edges)
    {
        for (dvec2 e : ej)
            out << e[0] << " " << e[1] << " ";
        out << std::endl;
    }
    out.close();

    std::vector<dvec2> shape;
    getTilingShape(shape, a_tiling, edges);

    std::vector<std::vector<dvec2>> polygons_v;

    Vector<T, 4> transf; TV2 xy;
    getTranslationUnitPolygon(polygons_v, shape, a_tiling, transf, unit, unit, xy);
    // getTranslationUnitPolygon(polygons_v, shape, a_tiling, transf, 8.0, 8.0, xy);
    // getTranslationUnitPolygon(polygons_v, shape, a_tiling, transf, 5.0, 5.0, xy);
    // getTranslationUnitPolygon(polygons_v, shape, a_tiling, transf, 4.0, 4.0, xy);


    csk::IsohedralTiling default_tiling( csk::tiling_types[ IH ] );

    int n_tiling_vtx = a_tiling.numVertices();
    Eigen::MatrixXd ipc_vtx(n_tiling_vtx, 2);
    Eigen::MatrixXd ipc_vtx_rest(n_tiling_vtx, 2);
    Eigen::MatrixXi ipc_edges(n_tiling_vtx, 2), ipc_faces;
    for (int i = 0; i < a_tiling.numVertices(); i++)
    {
        glm::vec2 vtx = a_tiling.getVertex(i);
        glm::vec2 vtx_rest = default_tiling.getVertex(i);
        ipc_vtx(i, 0) = vtx[0]; ipc_vtx(i, 1) = vtx[1];
        ipc_vtx_rest(i, 0) = vtx_rest[0]; ipc_vtx_rest(i, 1) = vtx_rest[1];
        ipc_edges(i, 0) = i; ipc_edges(i, 1) = (i+1)%n_tiling_vtx; 
    }

    auto signedArea = [&](const Eigen::MatrixXd& vtx) -> T
    {
        T area = 0.0;
        for (int i = 0; i < n_tiling_vtx; i++)
        {
            int vi = ipc_edges(i, 0), vj = ipc_edges(i, 1);
            TV xi = vtx.row(vi), xj = vtx.row(vj);
            area += 0.5 * (xi[0] * xj[1] - xi[1] * xj[0]);
        }
        
        return area;
    };

    T area_rest = std::abs(signedArea(ipc_vtx_rest));
    T area = std::abs(signedArea(ipc_vtx));
    // std::cout << std::abs(area/area_rest) * 100.0 << "%" << std::endl;
    if (std::abs(area/area_rest) < 0.15)
        return false;

    // if (ipc::has_intersections(ipc_vtx, ipc_edges, ipc_faces))
    //     return false;
    
    // std::ofstream obj_out("tiling_unit.obj");
    // for (int i = 0; i < a_tiling.numVertices(); i++)
    // {
    //     glm::vec2 vtx = a_tiling.getVertex(i);
    //     obj_out << "v " << vtx[0] << " " << vtx[1] << " 0.0" << std::endl;
    // }
    // for (int i = 0; i < a_tiling.numVertices(); i++)
    // {
    //     obj_out << "l " << i+1 << " " << (i+1)%a_tiling.numVertices() + 1 << std::endl;
    // }
    // obj_out.close();
    // std::exit(0);
    

    Vector<T, 8> periodic;
    periodic.head<2>() = TV2(0,0);
    periodic.segment<2>(2) = periodic.head<2>() + T(n_unit) * TV2(transf[0],transf[2]);
    periodic.segment<2>(4) = periodic.head<2>() + T(n_unit) * TV2(transf[0],transf[2]) + T(n_unit) * TV2(transf[1],transf[3]);
    periodic.segment<2>(6) = periodic.head<2>() + T(n_unit) * TV2(transf[1],transf[3]);

    // TM R = rotMat(angle);
    TM2 R = TM2::Identity();

    ClipperLib::Paths polygons(polygons_v.size());
    T mult = 1e12;
    for(int i=0; i<polygons_v.size(); ++i)
    {
        for(int j=0; j<polygons_v[i].size(); ++j)
        {
            TV2 curr = R * TV2(polygons_v[i][j][0]-xy[0], polygons_v[i][j][1]-xy[1]);
            polygons[i] << ClipperLib::IntPoint(curr[0]*mult, curr[1]*mult);
            // polygons[i] << ClipperLib::IntPoint((polygons_v[i][j][0]-xy[0])*mult, 
            //     (polygons_v[i][j][1]-xy[1])*mult);
            // std::cout << " " << polygons[i][j];
        }
        // break;
        // std::cout << std::endl;
    }
    // std::ofstream polygon_obj("polygon_obj.obj");
    // // for (auto polygon : polygons)
    // for (int i = 0; i < polygons.size(); i++)
    // {
    //     auto polygon = polygons[i];
        
    //     for (auto vtx : polygon)
    //         polygon_obj << "v " << vtx.X << " " << vtx.Y << " 0" << std::endl;
        
    // }
    // polygon_obj.close();

    periodic.segment<2>(2) = R * periodic.segment<2>(2);
    periodic.segment<2>(4) = R * periodic.segment<2>(4);
    periodic.segment<2>(6) = R * periodic.segment<2>(6);
    
    T distance = -1.5;
    ClipperLib::Paths final_shape;

    ClipperLib::ClipperOffset c;
    c.AddPaths(polygons, ClipperLib::jtSquare, ClipperLib::etClosedPolygon);
    
    c.Execute(final_shape, distance*mult);
    // saveClip(final_shape, periodic, mult, "tiling_unit_clip_in_x.obj", true);
    shapeToPolygon(final_shape, eigen_polygons, mult);
    periodicToBase(periodic, eigen_base);
    return true;
}

// ============================ PRIVATE ============================
glm::dmat3 IsohedralTiling::centrePSRect(T xmin, T ymin, T xmax, T ymax)
{
    T sc = std::min( 6.5*72.0 / (xmax-xmin), 9.0*72.0 / (ymax-ymin) );
    return glm::dmat3( 1, 0, 0, 0, 1, 0, 4.25*72.0, 5.5*72.0, 1.0 )
        * glm::dmat3( sc, 0, 0, 0, sc, 0, 0, 0, 1 )
        * glm::dmat3( 1, 0, 0, 0, 1, 0, -0.5*(xmin+xmax), -0.5*(ymin+ymax), 1 );
}

std::vector<glm::dvec2> IsohedralTiling::outShapeVec(const std::vector<glm::dvec2>& vec, const glm::dmat3& M)
{
    std::vector<glm::dvec2> data_points;

    glm::dvec2 p = M * glm::dvec3( vec.back(), 1.0 );
    data_points.push_back(glm::dvec2(p[0], p[1]));

    for( size_t idx = 0; idx < vec.size(); idx += 3 ) {
        glm::dvec2 p1 = M * glm::dvec3( vec[idx], 1.0 );
        glm::dvec2 p2 = M * glm::dvec3( vec[idx+1], 1.0 );
        glm::dvec2 p3 = M * glm::dvec3( vec[idx+2], 1.0 );

        data_points.push_back(glm::dvec2(p1[0], p1[1]));
        data_points.push_back(glm::dvec2(p2[0], p2[1]));
        data_points.push_back(glm::dvec2(p3[0], p3[1]));
    }

    return data_points;
}

void IsohedralTiling::getTilingShape(std::vector<dvec2>& shape, const csk::IsohedralTiling& tiling, 
    const std::vector<std::vector<dvec2>>& edges)
{
    for( auto i : tiling.shape() ) {
        // Get the relevant edge shape created above using i->getId().
        const std::vector<dvec2>& ed = edges[ i->getId() ];
        // Also get the transform that maps to the line joining consecutive
        // tiling vertices.
        const glm::dmat3& TT = i->getTransform();

        // If i->isReversed() is true, we need to run the parameterization
        // of the path backwards.
        if( i->isReversed() ) {
            for( size_t idx = 1; idx < ed.size(); ++idx ) {
                shape.push_back( TT * glm::dvec3( ed[ed.size()-1-idx], 1.0 ) );
            }
        } else {
            for( size_t idx = 1; idx < ed.size(); ++idx ) {
                shape.push_back( TT * glm::dvec3( ed[idx], 1.0 ) );
            }
        }
    }
}

void IsohedralTiling::getTilingEdges(const csk::IsohedralTiling& tiling,
        const Vector<T, 4>& eij,
        std::vector<std::vector<dvec2>>& edges)
{
    using namespace csk;
    for( U8 idx = 0; idx < tiling.numEdgeShapes(); ++idx ) {
        std::vector<dvec2> ej;

        ej.push_back( dvec2( 0, 0.0 ) );
        ej.push_back( dvec2( eij[0], eij[1] ) );
        ej.push_back( dvec2( eij[2], eij[3] ) );
        ej.push_back( dvec2( 1.0, 0.0 ) );
        
        // Now, depending on the edge shape class, enforce symmetry 
        // constraints on edges.
        switch( tiling.getEdgeShape( idx ) ) {
        case J: 
            break;
        case U:
            ej[2].x = 1.0 - ej[1].x;
            ej[2].y = ej[1].y;
            break;
        case S:
            ej[2].x = 1.0 - ej[1].x;
            ej[2].y = -ej[1].y;
            break;
        case I:
            ej[1].y = 0.0;
            ej[2].y = 0.0;
            break;
        }
        edges[idx] = ej;   
    }
}

void IsohedralTiling::shapeToPolygon(ClipperLib::Paths& final_shape, 
    std::vector<std::vector<TV2>>& polygons, T mult)
{
    polygons.resize(final_shape.size());

    for(int i=0; i<final_shape.size(); ++i)
    {
        for(int j=0; j<final_shape[i].size(); ++j)
        {
            TV2 cur_point = TV2(final_shape[i][j].X/mult, final_shape[i][j].Y/mult);

            if(j==final_shape[i].size()-1)
            {
                if((cur_point-polygons[i].front()).norm()>1e-4)
                    polygons[i].push_back(cur_point);
            }
            else if(j>0)
            {
                if((cur_point-polygons[i].back()).norm()>1e-4)
                    polygons[i].push_back(cur_point);
            }
            else
                polygons[i].push_back(cur_point);	
        }
    }
}

void IsohedralTiling::periodicToBase(const Vector<T, 8>& periodic, std::vector<TV2>& eigen_base)
{
    eigen_base.resize(0);

    eigen_base.push_back(TV2(periodic[0], periodic[1]));
    eigen_base.push_back(TV2(periodic[2], periodic[3]));
    eigen_base.push_back(TV2(periodic[4], periodic[5]));
    eigen_base.push_back(TV2(periodic[6], periodic[7]));
}

void IsohedralTiling::getTranslationUnitPolygon(std::vector<std::vector<dvec2>>& polygons_v,
        const std::vector<dvec2>& shape, const csk::IsohedralTiling& tiling, 
        Vector<T, 4>& transf, int width, int depth, TV2& xy_shift)
{
    int min_y=10000, max_y=-10000, min_x=10000, max_x=-10000;
    
    int ii=0;
    int extension = 4;
    dmat3 M = centrePSRect( -width, -depth, width, depth );
    for( auto i : tiling.fillRegion( -extension * width, -extension * depth, extension * width, extension * depth ) ) 
    {
        dmat3 T = M * i->getTransform();

        std::vector<dvec2> outv = outShapeVec( shape, T );

        if(T[0][0]!=T[1][1])
            std::reverse(outv.begin(), outv.end());

        min_y = std::min(i->getT2(), min_y);
        max_y = std::max(i->getT2(), max_y);

        min_x = std::min(i->getT1(), min_x);
        max_x = std::max(i->getT1(), max_x);

        polygons_v.push_back(outv);
    }

    int chosen_x = (max_x+min_x)/2;
    int chosen_y = (max_y+min_y)/2;

    TV2 xy, x1y, xy1, x1y1;

    for( auto i : tiling.fillRegion( -extension * width, -extension * depth, extension * width, extension * depth ) ) {

        dmat3 T = M * i->getTransform();
        std::vector<dvec2> outv = outShapeVec( shape, T );
    
        if(T[0][0]!=T[1][1])
            std::reverse(outv.begin(), outv.end());

        if(i->getT1() == chosen_x && i->getT2() == chosen_y && i->getAspect()==0)
            xy << outv[0].x, outv[0].y;		

        if(i->getT1() == chosen_x+1 && i->getT2() == chosen_y && i->getAspect()==0)
            x1y << outv[0].x, outv[0].y;		

        if(i->getT1() == chosen_x && i->getT2() == chosen_y+1 && i->getAspect()==0)
            xy1 << outv[0].x, outv[0].y;		

        if(i->getT1() == chosen_x+1 && i->getT2() == chosen_y+1 && i->getAspect()==0)
            x1y1 << outv[0].x, outv[0].y;		

    }

    transf.setZero();
    transf.head<2>() = x1y - xy;
    transf.tail<2>() = xy1 - xy;
    
    T temp1 = transf[0];
    T temp2 = transf[3];

    transf[0] = transf[2];
    transf[3] = transf[1];

    transf[1] = temp2;
    transf[2] = temp1;

    if(transf[0]<0)
        transf.head<2>() *= -1;
    if(transf[3]<0)
        transf.tail<2>() *= -1;

    T temp = transf[2];
    transf[2] = transf[1];
    transf[1] = temp;
    xy_shift = xy;
}

void IsohedralTiling::saveClip(const ClipperLib::Paths& final_shape, 
        const Vector<T, 8>& periodic, T mult,
        const std::string& filename, bool add_box)
{
    std::ofstream out(filename);
    if (add_box)
        for (int i = 0; i < 4; i++)
            out << "v " <<  periodic.segment<2>(i * 2).transpose() << " 0" << std::endl;
    for (auto polygon : final_shape)
    {
        for (auto vtx : polygon)
        {
            out << "v " << vtx.X / mult << " " << vtx.Y / mult << " 0" << std::endl;
        }
    }
    int cnt = 5;
    if (add_box)
    {
        out << "l 1 2" << std::endl;
        out << "l 2 3" << std::endl;
        out << "l 3 4" << std::endl;
        out << "l 4 1" << std::endl;
    }
    else
        cnt = 1;
    for (auto polygon : final_shape)
    {
        for (int i = 0; i < polygon.size(); i++)
        {
            int j = (i + 1) % polygon.size();
            out << "l " << cnt + i << " " << cnt + j << std::endl;
        }
        cnt += polygon.size();
    }
    out.close();
}