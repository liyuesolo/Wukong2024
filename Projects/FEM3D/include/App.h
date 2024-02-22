#ifndef APP_H
#define APP_H

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

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

public:
    void initializeScene()
    {
        polyscope::options::autocenterStructures = true;
        polyscope::view::windowWidth = 3000;
        polyscope::view::windowHeight = 2000;
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
        polyscope::options::groundPlaneHeightFactor = 0.; 
        polyscope::options::shadowDarkness = 0.4;

        // Initialize polyscope
        polyscope::init();
        
        psMesh = polyscope::registerSurfaceMesh("surface mesh", simulation.surface_vertices, simulation.surface_indices);
        psMesh->setSmoothShade(false);
        psMesh->setSurfaceColor(glm::vec3(0.255, 0.514, 0.996));
        psMesh->setEdgeWidth(1.0);

        polyscope::state::userCallback = [&](){ sceneCallback(); };
    }

    void sceneCallback()
    {
        if (ImGui::Button("RunSim")) 
        {
            run_sim = true;
        }
        if (ImGui::Button("Animate Modes")) 
        {
            animate_modes = true;
        }
        if (ImGui::Button("Next Modes")) 
        {
            modes = (modes + 1) % eigen_values.rows();
        }
        if (ImGui::Button("stop")) 
        {
            animate_modes = false;
            run_sim = false;
        }
        if (ImGui::Button("Compute Linear Modes")) 
        {
            simulation.computeLinearModes(eigen_vectors, eigen_values);   
        }
        if (ImGui::Button("staticSolveOneStep")) 
        {
            simulation.advanceOneStep(static_solve_step++);
            psMesh = polyscope::registerSurfaceMesh("surface mesh", simulation.surface_vertices, simulation.surface_indices);
        }
        if (animate_modes && !run_sim)
        {
            t += 0.1;
            simulation.deformed = simulation.undeformed + simulation.u + eigen_vectors.col(modes) * std::sin(t);
            simulation.updateSurfaceVertices();
            
            psMesh->updateVertexPositions(simulation.surface_vertices);
        }
        if (!animate_modes && run_sim)
        {
            bool finished = simulation.advanceOneStep(static_solve_step++);
            psMesh->updateVertexPositions(simulation.surface_vertices);
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