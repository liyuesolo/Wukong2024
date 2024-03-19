#ifndef APP_H
#define APP_H

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "VecMatDef.h"
#include "Util.h"
#include "DiscreteShellHomogenization.h"
#include "Homogenization.h"

enum VisualizationUnit
{
    Triangle, Rod, ImplicitSurface
};

template<class Simulation>
class App
{
public:
    Simulation& simulation;
    int static_solve_step = 0;

    polyscope::SurfaceMesh* psMesh;

    T t = 0;

    bool animate_modes = false;
    bool run_sim = false;
    int modes = 0;
    Eigen::MatrixXd eigen_vectors;
    Eigen::VectorXd eigen_values;

    Eigen::MatrixXd meshV;
    Eigen::MatrixXi meshF;

    VisualizationUnit visualization_unit = Triangle;
public:
    void initializeScene();

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
        
        if (ImGui::Button("staticSolveOneStep")) 
        {
            // simulation.advanceOneStep(static_solve_step++);
            // Eigen::MatrixXd meshV;
            // Eigen::MatrixXi meshF;
            // vectorToIGLMatrix<T, 3>(simulation.deformed, meshV);
            // vectorToIGLMatrix<int, 3>(simulation.faces, meshF);
            // polyscope::registerSurfaceMesh("surface mesh", meshV, meshF);
        }
        if (!animate_modes && run_sim)
        {
            // bool finished = simulation.advanceOneStep(static_solve_step++);
            // Eigen::MatrixXd meshV;
            // vectorToIGLMatrix<T, 3>(simulation.deformed, meshV);
            // psMesh->updateVertexPositions(meshV);
            // if (finished)
            //     run_sim = false;
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