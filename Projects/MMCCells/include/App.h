#ifndef APP_H
#define APP_H

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/implicit_helpers.h"

template<class Simulation>
class App
{
public:
    using TV = Vector<T, 3>;
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

    bool show_levelset = true;

public:
    void initializeScene()
    {
        polyscope::options::autocenterStructures = true;
        polyscope::view::windowWidth = 3000;
        polyscope::view::windowHeight = 2000;
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
        polyscope::options::groundPlaneHeightFactor = 0.; 
        polyscope::options::shadowDarkness = 0.4;

        // Initialize polyscope
        polyscope::init();
        if (show_levelset)
        {
            polyscope::ImplicitRenderOpts opts;
            opts.subsampleFactor = 20;
            polyscope::ImplicitRenderMode mode = polyscope::ImplicitRenderMode::SphereMarch;
            // for (int i = 0; i < simulation.n_components; i++)
            // {
            //     auto sdf_i = [&](glm::vec3 p) { return simulation.levelsetValue(TV(p[0], p[1], p[2]), i); };
            //     polyscope::DepthRenderImageQuantity* img = polyscope::renderImplicitSurface("sdf_i", sdf_i, mode, opts);
            // }
            auto sdf_0 = [&](glm::vec3 p) { return simulation.levelsetValue(TV(p[0], p[1], p[2])); };
            polyscope::DepthRenderImageQuantity* img0 = polyscope::renderImplicitSurface("sdf0", sdf_0, mode, opts);
            
        }

        // psMesh = polyscope::registerSurfaceMesh("surface mesh", simulation.surface_vertices, simulation.surface_indices);
        // psMesh->setSmoothShade(false);
        // psMesh->setSurfaceColor(glm::vec3(0.255, 0.514, 0.996));
        // psMesh->setEdgeWidth(1.0);

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
        if (ImGui::Button("show levelset")) 
        {
            show_levelset = !show_levelset;
        }
        if (ImGui::Button("staticSolveOneStep")) 
        {
            // simulation.advanceOneStep(static_solve_step++);
            // psMesh = polyscope::registerSurfaceMesh("surface mesh", simulation.surface_vertices, simulation.surface_indices);
        }
        if (ImGui::Button("DerivativeTest")) 
        {
            // simulation.checkTotalGradientScale(true);
        }
        if (animate_modes && !run_sim)
        {
            t += 0.1;
            // simulation.deformed = simulation.undeformed + simulation.u + eigen_vectors.col(modes) * std::sin(t);
            // simulation.updateSurfaceVertices();
            
            // psMesh->updateVertexPositions(simulation.surface_vertices);
        }
        if (!animate_modes && run_sim)
        {
            // bool finished = simulation.advanceOneStep(static_solve_step++);
            // psMesh->updateVertexPositions(simulation.surface_vertices);
            // if (finished)
            //     run_sim = false;
        }
        if (show_levelset)
        {
            polyscope::ImplicitRenderOpts opts;
            opts.subsampleFactor = 20;
            polyscope::ImplicitRenderMode mode = polyscope::ImplicitRenderMode::SphereMarch;
            auto sdf_0 = [&](glm::vec3 p) { return simulation.levelsetValue(TV(p[0], p[1], p[2])); };
            polyscope::DepthRenderImageQuantity* img0 = polyscope::renderImplicitSurface("sdf0", sdf_0, mode, opts);
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