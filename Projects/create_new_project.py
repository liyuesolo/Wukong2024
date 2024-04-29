import sys
import os

def createCMakeLists(project_name, use_Spectra = True, use_IPC = False):
    cmake_file = project_name + "/CMakeLists.txt"
    if not os.path.exists(cmake_file):
        f = open(cmake_file, "w+")
        f.write("set(DEPS boost_filesystem tbb \n \
            tbbmalloc tbbmalloc_proxy \n \
            mkl_intel_lp64 igl::core \n \
            mkl_sequential mkl_core  \n \
            suitesparseconfig metis cholmod amd camd ccolamd colamd spqr \n \
            gmp mpfr)\n\n"
            )
        f.write("file(GLOB HEADERS \"include/*.h\" \"autodiff/*.h\")\n")
        f.write("file(GLOB SOURCES \"src/*.cpp\" \"autodiff/*.cpp\")\n\n")

        f.write("add_executable("+project_name+" ${HEADERS} ${SOURCES})\n")
        if use_Spectra:
            f.write("target_include_directories("+project_name+" PUBLIC ../../Libs/spectra/include)\n")    
        f.write("target_include_directories("+project_name+" PUBLIC ../../Deps/libigl/include)\n")
        f.write("target_link_libraries("+project_name+" ${DEPS} polyscope)")
        f.close()

def createTemplateFiles(base_folder, project_name):
    os.system("cp DiscreteShell/include/VecMatDef.h " + base_folder + "/include")
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
        f.write("\t" + project_name + "() {} \n")
        f.write("\t~" + project_name + "() {} \n")
        f.write("};\n")
        f.write("\n\n")
        f.write("#endif\n")
        f.close()
    cpp_file = base_folder + "/src/" + project_name + ".cpp"
    if not os.path.exists(cpp_file):
        f = open(cpp_file, "w+")
        f.write("#include \"../include/" + project_name + ".h\"")
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