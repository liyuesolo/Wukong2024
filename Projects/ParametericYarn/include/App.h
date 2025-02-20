#ifndef APP_H 
 #define APP_H 

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

template<class Simulation>  
class App  
{  
public:  
    Simulation& simulation;  
    int static_solve_step = 0;

    // polyscope::SurfaceMesh* mesh;
    
    polyscope::PointCloud* point_cloud;

    T t = 0;

    bool animate_modes = false;
    bool run_sim = false;
    int modes = 0;
    Eigen::MatrixXd eigen_vectors;
    Eigen::VectorXd eigen_values;

    Eigen::MatrixXd meshV;
    Eigen::MatrixXi meshF;
public:
    void initializeScene()
    {
        polyscope::options::programName = "WuKongSim";
        polyscope::options::autocenterStructures = true;
        polyscope::view::windowWidth = 3000;
        polyscope::view::windowHeight = 2000;
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
        polyscope::options::groundPlaneHeightFactor = 0.6; 
        polyscope::options::shadowDarkness = 0.4;
        // Initialize polyscope
        polyscope::init();
        vectorToIGLMatrix<T, 3>(simulation.deformed, meshV);
        //vectorToIGLMatrix<int, 3>(simulation.faces, meshF);
        point_cloud = polyscope::registerPointCloud("points", meshV);
        //psMesh->setSmoothShade(false);
        point_cloud->setPointColor(glm::vec3(0.255, 0.514, 0.996));
        //psMesh->setEdgeWidth(1.0);        
        polyscope::state::userCallback = [&](){ sceneCallback(); };
    }
    void sceneCallback() {}
        void run()
    {
        polyscope::show();
    }

public:
    App(Simulation& sim) : simulation(sim) {}
    ~App() {}
};

#endif
