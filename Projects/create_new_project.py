import sys
import os

def createCMakeLists(project_name, use_Spectra = True, use_IPC = False):
    cmake_file = project_name + "/CMakeLists.txt"
    if not os.path.exists(cmake_file):
        f = open(cmake_file, "w+")
        f.write("if(APPLE)\n")
        f.write("\tset(DEPS igl::core)\n")
        f.write("else()\n")
        f.write("\tset(DEPS boost_filesystem tbb \n \
            tbbmalloc tbbmalloc_proxy \n \
            mkl_intel_lp64 igl::core \n \
            mkl_sequential mkl_core  \n \
            suitesparseconfig metis cholmod amd camd ccolamd colamd spqr \n \
            gmp mpfr)\n\n"
            )
        f.write("endif()\n\n")
        f.write("file(GLOB HEADERS \"include/*.h\" \"autodiff/*.h\")\n")
        f.write("file(GLOB SOURCES \"src/*.cpp\" \"autodiff/*.cpp\")\n\n")

        f.write("add_executable("+project_name+" ${HEADERS} ${SOURCES})\n")
        f.write("if(APPLE)\n")
        f.write(" \
    find_package(CHOLMOD REQUIRED) \n \
    include_directories(${CHOLMOD_INCLUDES}) \n \
    find_package(TBB REQUIRED) \n")
        f.write("\ttarget_link_libraries("+project_name+" ${CHOLMOD_LIBRARIES} TBB::tbb)\n")
        f.write("endif()\n")
        if use_Spectra:
            f.write("target_include_directories("+project_name+" PUBLIC ../../Libs/spectra/include)\n")    
        f.write("target_include_directories("+project_name+" PUBLIC ../../Deps/libigl/include)\n")
        f.write("target_link_libraries("+project_name+" ${DEPS} polyscope)")
        f.close()

def createTemplateFiles(base_folder, project_name):
    os.system("cp DiscreteShell/include/VecMatDef.h " + base_folder + "/include")
    os.system("cp DiscreteShell/include/Timer.h " + base_folder + "/include")
    os.system("cp DiscreteShell/include/Util.h " + base_folder + "/include")
    if not os.path.exists(base_folder + "/include/" + project_name + ".h"):
        f = open(base_folder + "/include/" + project_name + ".h", "w+")
        f.write("#ifndef " + project_name.upper() + "_H\n")
        f.write("#define " + project_name.upper() + "_H\n")
        f.write("\n\n")
        f.write("#include <utility> \n\
#include <iostream>\n\
#include <fstream>\n\
#include <Eigen/Geometry>\n\
#include <Eigen/Core>\n\
#include <Eigen/Sparse>\n\
#include <Eigen/Dense>\n\
#include <tbb/tbb.h>\n\
#include <unordered_map>\n\
#include <iomanip>\n\n\n")
        f.write("#include \"VecMatDef.h\"\n\n")
        f.write("#include \"Timer.h\"\n\n")

        f.write("class " + project_name + "\n")
        f.write("{\npublic:\n")
        f.write("\tusing VectorXT = Matrix<T, Eigen::Dynamic, 1>;\n")
        f.write("\tusing MatrixXT = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;\n")
        f.write("\tusing VectorXi = Vector<int, Eigen::Dynamic>;\n")
        f.write("\tusing MatrixXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;\n")
        f.write("\tusing VtxList = std::vector<int>;\n")
        f.write("\tusing StiffnessMatrix = Eigen::SparseMatrix<T>;\n")
        f.write("\tusing Entry = Eigen::Triplet<T>;\n")
        f.write("\tusing TV = Vector<T, 3>;\n")
        f.write("\tusing TV2 = Vector<T, 2>;\n")
        f.write("\tusing TM2 = Matrix<T, 2, 2>;\n")
        f.write("\tusing TV3 = Vector<T, 3>;\n")
        f.write("\tusing IV = Vector<int, 3>;\n")
        f.write("\tusing IV2 = Vector<int, 2>;\n")
        f.write("\tusing TM = Matrix<T, 3, 3>;\n")
        f.write("\n")
        f.write("public:\n")
        f.write("\tbool verbose = false;\n")
        f.write("\tbool dynamics = true;\n")
        f.write("\tVectorXT deformed, undeformed, u;\n")
        f.write("\tVectorXT external_force;\n")
        f.write("\tbool run_diff_test = false;\n")
        f.write("\tstd::unordered_map<int, T> dirichlet_data;\n")
        f.write("\tT newton_tol = 1e-6;\n")
        f.write("\tint max_newton_iter = 100;\n")
        f.write("\n")
        f.write("public:\n")
        f.write(" \
    template <class OP> \n \
    void iterateDirichletDoF(const OP& f) \n \
    { \n \
        for (auto dirichlet : dirichlet_data) \n \
        { \n \
            f(dirichlet.first, dirichlet.second); \n \
        } \n \
    }\n")
        f.write("\n")
        f.write("public:\n")
        f.write("\t" + project_name + "() {} \n")
        f.write("\t~" + project_name + "() {} \n")
        f.write("\n")
        f.write("\tvoid initialize() {}\n")
        f.write("\tbool advanceOneStep(int step);\n")
        f.write("\tbool advanceOneTimeStep();\n")
        f.write("\tbool linearSolve(StiffnessMatrix& K, const VectorXT& residual,VectorXT& du);\n")
        f.write("\tvoid buildSystemMatrix(StiffnessMatrix& K);\n")
        f.write("\tT lineSearchNewton(const VectorXT& residual);\n")
        f.write("\tT computeTotalEnergy();\n")
        f.write("\tT computeResidual(VectorXT& residual);\n")
        f.write("\tvoid projectDirichletDoFMatrix(StiffnessMatrix& A,\
                                    const std::unordered_map<int, T>& data);\n")
        f.write("};\n")
        f.write("\n\n")
        f.write("#endif\n")
        f.close()
    cpp_file = base_folder + "/src/" + project_name + ".cpp"
    if not os.path.exists(cpp_file):
        f = open(cpp_file, "w+")
        f.write("#include <Eigen/CholmodSupport>\n")
        f.write("#include \"../include/" + project_name + ".h\"")
        f.write("\n\n")
        f.write("T " + project_name + "::computeTotalEnergy()\n")
        f.write("{\n")
        f.write(" \n \
    if (!run_diff_test) \n \
        iterateDirichletDoF([&](int offset, T target) \n \
                            { u[offset] = target; }); \n \
    deformed = undeformed + u; \n \
    T energy = 0.0; \n \
    \n \
    // if (dynamics) \n \
    //    addInertialEnergy(energy); \n \
    //energy -= u.dot(external_force); \n \
    return energy;\n \
        ")
        f.write("}\n")
        f.write("\n")
        f.write("T " + project_name + "::computeResidual(VectorXT& residual)\n")
        f.write("{ \n \
    if (!run_diff_test) \n \
        iterateDirichletDoF([&](int offset, T target) \n \
                            { u[offset] = target; }); \n \
    deformed = undeformed + u; \n \
    \n \
    if (!run_diff_test) \n \
        iterateDirichletDoF([&](int offset, T target) \n \
                            { residual[offset] = 0; }); \n \
    return residual.norm(); \n")
        f.write("}\n")
        f.write("\n")
        f.write("bool " + project_name + "::linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du)\n")
        f.write("{\n")
        f.write(" \n \
    START_TIMING(LinearSolve) \n \
    Eigen::CholmodSupernodalLLT<StiffnessMatrix, Eigen::Lower> solver; \n \
 \n \
    T alpha = 1e-6; \n \
    if (!dynamics) \n \
    { \n \
        StiffnessMatrix H(K.rows(), K.cols()); \n \
        H.setIdentity(); \n \
        H.diagonal().array() = 1e-10; \n \
        K += H; \n \
    } \n \
    solver.analyzePattern(K); \n \
 \n \
    int indefinite_count_reg_cnt = 0, invalid_search_dir_cnt = 0, \n \
        invalid_residual_cnt = 0; \n \
    int i = 0; \n \
    T dot_dx_g = 0.0; \n \
    for (; i < 50; i++) \n \
    { \n \
        solver.factorize(K); \n \
        if (solver.info() == Eigen::NumericalIssue) \n \
        { \n \
            K.diagonal().array() += alpha; \n \
            alpha *= 10; \n \
            indefinite_count_reg_cnt++; \n \
            continue; \n \
        } \n \
 \n \
        du = solver.solve(residual); \n \
 \n \
        dot_dx_g = du.normalized().dot(residual.normalized()); \n \
 \n \
        int num_negative_eigen_values = 0; \n \
        int num_zero_eigen_value = 0; \n \
 \n \
        bool positive_definte = num_negative_eigen_values == 0; \n \
        bool search_dir_correct_sign = dot_dx_g > 1e-6; \n \
        if (!search_dir_correct_sign) \n \
        { \n \
            invalid_search_dir_cnt++; \n \
        } \n \
 \n \
      \n \
        bool solve_success = du.norm() < 1e3; \n \
 \n \
        if (!solve_success) \n \
            invalid_residual_cnt++; \n \
         \n \
 \n \
        if (positive_definte && search_dir_correct_sign && solve_success) \n \
        { \n \
 \n \
            if (verbose) \n \
            { \n \
                std::cout << \"\t===== Linear Solve ===== \" << std::endl; \n \
                std::cout << \"\tnnz: \" << K.nonZeros() << std::endl; \n \
                 \n \
                std::cout << \"\t# regularization step \" << i << \" indefinite \" \n \
                          << indefinite_count_reg_cnt << \" invalid search dir \" \n \
                          << invalid_search_dir_cnt << \" invalid solve \" \n \
                          << invalid_residual_cnt << std::endl; \n \
                std::cout << \"\tdot(search, -gradient) \" << dot_dx_g \n \
                          << std::endl; \n \
                std::cout << \"\t======================== \" << std::endl; \n \
                FINISH_TIMING_PRINT(LinearSolve) \n \
            } \n \
            return true; \n \
        } \n \
        else \n \
        { \n \
            K.diagonal().array() += alpha; \n \
            alpha *= 10; \n \
        } \n \
    } \n \
    if (verbose) \n \
    { \n \
        std::cout << \"\t===== Linear Solve ===== \" << std::endl; \n \
        std::cout << \"\tnnz: \" << K.nonZeros() << std::endl; \n \
        // std::cout << \"\t takes \" << t.elapsed_sec() << \"s\" << std::endl; \n \
        std::cout << \"\t# regularization step \" << i << \" indefinite \" \n \
                  << indefinite_count_reg_cnt << \" invalid search dir \" \n \
                  << invalid_search_dir_cnt << \" invalid solve \" \n \
                  << invalid_residual_cnt << std::endl; \n \
        std::cout << \"\tdot(search, -gradient) \" << dot_dx_g << std::endl; \n \
        std::cout << \"\t======================== \" << std::endl; \n \
        FINISH_TIMING_PRINT(LinearSolve) \n \
    } \n \
        return false;\n ")
        f.write("}\n")
        f.write("\n")
        f.write("void " + project_name + "::buildSystemMatrix(StiffnessMatrix& K)\n")
        f.write("{\n")
        f.write(" \
    int n_dof = deformed.rows(); \n \
    if (!run_diff_test) \n \
        iterateDirichletDoF([&](int offset, T target) \n \
                            { u[offset] = target; }); \n \
    deformed = undeformed + u; \n \
    std::vector<Entry> entries; \n \
    \n \
    K.resize(n_dof, n_dof); \n \
    K.setFromTriplets(entries.begin(), entries.end()); \n \
\n \
    if (!run_diff_test) \n \
        projectDirichletDoFMatrix(K, dirichlet_data); \n \
    K.makeCompressed();\n ")
        f.write("}\n")
        f.write("\n")
        f.write("T " + project_name + "::lineSearchNewton(const VectorXT& residual)\n")
        f.write("{\n")
        f.write(" \
    VectorXT du = residual; \n \
    du.setZero(); \n \
\n \
    du = residual; \n \
    StiffnessMatrix K(residual.rows(), residual.rows()); \n \
    buildSystemMatrix(K); \n \
    bool success = linearSolve(K, residual, du); \n \
    if (!success) \n \
    {   \n \
        std::cout << \"Linear Solve Failed\" << std::endl; \n \
        return 1e16; \n \
    } \n \
\n \
    T norm = du.norm(); \n \
    if (verbose) \n \
        std::cout << \"\t|du | \" << norm << std::endl; \n \
\n \
    T E0 = computeTotalEnergy(); \n \
    T alpha = 1.0; \n \
    int cnt = 0; \n \
    VectorXT u_current = u; \n \
    while (true) \n \
    { \n \
        u = u_current + alpha * du; \n \
        T E1 = computeTotalEnergy(); \n \
        if (E1 - E0 < 0 || cnt > 10) \n \
        { \n \
            break; \n \
        } \n \
        alpha *= 0.5; \n \
        cnt += 1; \n \
    } \n \
    return alpha * du.norm(); \n")
        f.write("}\n")
        f.write("\n")
        f.write("bool " + project_name + "::advanceOneStep(int step)\n")
        f.write("{\n")
        f.write("\treturn false;\n")
        f.write("}\n")
        f.write("\n")
        f.write("bool " + project_name + "::advanceOneTimeStep()\n")
        f.write("{\n")
        f.write(" \n \
    int iter = 0;\n \
    while (true)\n \
    {\n \
        VectorXT residual = external_force;\n \
        T residual_norm = computeResidual(residual);\n \
        // residual_norms.push_back(residual_norm);\n \
        if (verbose)\n \
            std::cout << \"[NEWTON] iter \" << iter << \"/\" << max_newton_iter\n \
                      << \": residual_norm \" << residual_norm\n \
                      << \" tol: \" << newton_tol << std::endl;\n \
        if (residual_norm < newton_tol || iter == max_newton_iter)\n \
        {\n \
            std::cout << \"[NEWTON] iter \" << iter << \"/\" << max_newton_iter\n \
                      << \": residual_norm \" << residual_norm\n \
                      << \" tol: \" << newton_tol << std::endl;\n \
            break;\n \
        }\n \
        T du_norm = 1e10;\n \
        du_norm = lineSearchNewton(residual);\n \
        if (du_norm < 1e-10)\n \
            break;\n \
        iter++;\n \
    }\
    return true;\n ")
        f.write("}\n")
        f.write("void " + project_name + "::projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data)\n")
        f.write("{\n")
        f.write(" \
        for (auto iter : data) \n \
    { \n \
        A.row(iter.first) *= 0.0; \n \
        A.col(iter.first) *= 0.0; \n \
        A.coeffRef(iter.first, iter.first) = 1.0; \n \
    } \n \
        ")
        f.write("}\n")  
        f.close()
    
    main_cpp =  base_folder + "/src/main.cpp"
    if not os.path.exists(main_cpp):
        f = open(main_cpp, "w+")
        f.write("#include \"../include/"+project_name+".h\"\n")
        f.write("#include \"../include/App.h\"\n")
        f.write("\nint main()\n{\n\n")
        f.write("\t" + project_name + " " + project_name.lower() + ";\n")
        f.write("\tApp<"+project_name+"> app(" + project_name.lower() + ");\n")

        f.write("\tapp.initializeScene();\n\tapp.run();\n\n\treturn 0;\n}")
        f.close()

    app_header = base_folder + "/include/App.h"
    if not os.path.exists(app_header):
        f = open(app_header, "w+")
        f.write("#ifndef APP_H \n \
#define APP_H \n\
\n\
#include \"polyscope/polyscope.h\"\n\
#include \"polyscope/surface_mesh.h\"\n\
\n\
template<class Simulation>  \n\
class App  \n\
{  \n\
public:  \n\
    Simulation& simulation;  \n\
    int static_solve_step = 0;\n\
\n\
    polyscope::SurfaceMesh* psMesh;\n\
\n\
    T t = 0;\n\
\n\
    bool animate_modes = false;\n\
    bool run_sim = false;\n\
    int modes = 0;\n\
    Eigen::MatrixXd eigen_vectors;\n\
    Eigen::VectorXd eigen_values;\n\
\n\
    Eigen::MatrixXd meshV;\n\
    Eigen::MatrixXi meshF;\n\
public:\n\
    void initializeScene()\n\
    {\n\
        polyscope::options::autocenterStructures = true;\n\
        polyscope::view::windowWidth = 3000;\n\
        polyscope::view::windowHeight = 2000;\n\
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;\n\
        polyscope::options::groundPlaneHeightFactor = 0.6; \n\
        polyscope::options::shadowDarkness = 0.4;\n\
        // Initialize polyscope\n\
        polyscope::init();\n\
        //vectorToIGLMatrix<T, 3>(simulation.deformed, meshV);\n\
        //vectorToIGLMatrix<int, 3>(simulation.faces, meshF);\n\
        //psMesh = polyscope::registerSurfaceMesh(\"surface mesh\", meshV, meshF);\n\
        //psMesh->setSmoothShade(false);\n\
        //psMesh->setSurfaceColor(glm::vec3(0.255, 0.514, 0.996));\n\
        //psMesh->setEdgeWidth(1.0);        \n\
        polyscope::state::userCallback = [&](){ sceneCallback(); };\n\
    }\n\
    void sceneCallback() {}\n\
        void run()\n\
    {\n\
        polyscope::show();\n\
    }\n\
\n\
public:\n\
    App(Simulation& sim) : simulation(sim) {}\n\
    ~App() {}\n\
};\n\
\n\
#endif\n")
        f.close()

def createSubFolders(base_folder):
    sub_folders = ["include", "autodiff", "src"]
    for folder in sub_folders:
        if not os.path.exists(base_folder + "/" + folder):
            os.mkdir(base_folder + "/" + folder)
    
def addToCMake(project_name):
    f = open("CMakeLists.txt", "r+")
    found = False
    for line in f.readlines():
        if project_name in line.strip():
            found = True
            break
    if not found:
        f.write("\nadd_subdirectory(" + project_name + ")")
    f.close()

def main(project_name):
    if not os.path.exists(project_name):
        os.mkdir(project_name)
    createSubFolders(project_name)
    createCMakeLists(project_name, use_Spectra=True, use_IPC=False)
    createTemplateFiles(project_name, project_name)
    addToCMake(project_name)

if __name__ == '__main__':
    project_name = sys.argv[1] if len(sys.argv) >= 2 else ''
    main(project_name)