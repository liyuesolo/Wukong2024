#ifndef APP_H
#define APP_H

#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"

template<class Simulation>
class App
{
public:
    using TV = Vector<T, 3>;
    Simulation& simulation;
    int static_solve_step = 0;

    polyscope::CurveNetwork* rod_network;

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
        polyscope::options::groundPlaneHeightFactor = 0.2; 
        polyscope::options::shadowDarkness = 0.4;

        // Initialize polyscope
        polyscope::init();
        // std::cout << simulation.Rods.size() << std::endl;
        std::vector<glm::vec3> nodes;
        std::vector<std::array<size_t, 2>> edges;
        int global_cnt = 0;
        for(auto& rod : simulation.Rods)
        {
            for (int idx : rod->indices)
            {
                TV node_i;
                rod->x(idx, node_i);
                nodes.push_back(glm::vec3(node_i[0], node_i[1], node_i[2]));
            }
            for (int i = 0; i < rod->indices.size()-1; i++)
            {
                edges.push_back({global_cnt + i, global_cnt + i + 1});
            }
            global_cnt += rod->indices.size();
        }
        rod_network = polyscope::registerCurveNetwork("rod network", nodes, edges);
        // psMesh = polyscope::registerSurfaceMesh("surface mesh", simulation.surface_vertices, simulation.surface_indices);
        // psMesh->setSmoothShade(false);
        rod_network->setColor(glm::vec3(0.255, 0.514, 0.996));
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
        if (ImGui::Button("staticSolveOneStep")) 
        {
            // simulation.advanceOneStep(static_solve_step++);
            // psMesh = polyscope::registerSurfaceMesh("surface mesh", simulation.surface_vertices, simulation.surface_indices);
        }
        if (ImGui::Button("DerivativeTest")) 
        {
            // simulation.checkTotalGradientScale(true);
        }
        if (!animate_modes && run_sim)
        {
            bool finished = simulation.advanceOneStep(static_solve_step++);
            std::vector<glm::vec3> nodes;
            for(auto& rod : simulation.Rods)
            {
                for (int idx : rod->indices)
                {
                    TV node_i;
                    rod->x(idx, node_i);
                    nodes.emplace_back(node_i[0], node_i[1], node_i[2]);
                }
            }
            rod_network->updateNodePositions(nodes);
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