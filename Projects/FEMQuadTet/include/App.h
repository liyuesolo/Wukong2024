#ifndef APP_H 
 #define APP_H 

#include "polyscope/polyscope.h"
#include "polyscope/volume_mesh.h"
#include "Util.h"

template<class Simulation>  
class App  
{  
public:  
    Simulation& simulation;  
    int static_solve_step = 0;

    polyscope::VolumeMesh* psMesh;

    T t = 0;

    bool animate_modes = false;
    bool run_sim = false;
    int modes = 0;
    Eigen::MatrixXd eigen_vectors;
    Eigen::VectorXd eigen_values;

    Eigen::MatrixXd meshV;
    
public:
    void initializeScene()
    {
        polyscope::options::autocenterStructures = true;
        polyscope::view::windowWidth = 3000;
        polyscope::view::windowHeight = 2000;
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
        // polyscope::options::groundPlaneHeightFactor = 0.6; 
        // polyscope::options::shadowDarkness = 0.4;
        // Initialize polyscope
        polyscope::init();
        vectorToIGLMatrix<T, 3>(simulation.deformed.segment(0, simulation.n_linear_nodes * 3), meshV);
        //vectorToIGLMatrix<int, 3>(simulation.faces, meshF);
        psMesh = polyscope::registerTetMesh("volume mesh", meshV, simulation.linear_tet_indices);
        psMesh->setColor(glm::vec3(0.255, 0.514, 0.996));
        psMesh->setEdgeWidth(1.5);        
        polyscope::state::userCallback = [&](){ sceneCallback(); };
    }
    void sceneCallback() 
    {
        if (ImGui::Button("RunSim")) 
        {
            run_sim = true;
        }
        if (ImGui::Button("stop")) 
        {
            animate_modes = false;
            run_sim = false;
        }
        if (!animate_modes && run_sim)
        {
            bool finished = simulation.advanceOneStep(static_solve_step++);
            vectorToIGLMatrix<T, 3>(simulation.deformed.segment(0, simulation.n_linear_nodes * 3), meshV);
            psMesh->updateVertexPositions(meshV);
            if (finished)
                run_sim = false;
        }
    }
        void run()
    {
        polyscope::show();
    }

public:
    App(Simulation& sim) : simulation(sim) {}
    ~App() {}
};

#endif
