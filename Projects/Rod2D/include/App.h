#ifndef APP_H
#define APP_H

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"

template <class Simulation>
class App
{
public:
    using TV = Vector<T, 2>;

    Simulation& simulation;
    int static_solve_step = 0;

    polyscope::CurveNetwork* rod;

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
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
        // polyscope::options::groundPlaneHeightFactor = 0.6;
        // polyscope::options::shadowDarkness = 0.4;
        // Initialize polyscope
        polyscope::init();
        std::vector<glm::vec3> nodes;
        std::vector<std::array<size_t, 2>> edges;
        for (int i = 0; i < simulation.deformed.rows() / 2; i++)
        {
            TV node_i = simulation.deformed.segment(i * 2, 2);
            nodes.emplace_back(node_i[0], node_i[1], 0);
        }
        
        for (int i = 0; i < simulation.rod_indices.size() / 2; i++)
        {
            edges.push_back({size_t(simulation.rod_indices[i * 2]),
                             size_t(simulation.rod_indices[i * 2 + 1])});
        }
        rod = polyscope::registerCurveNetwork("rod", nodes, edges);
        rod->setColor(glm::vec3(0.255, 0.514, 0.996));
        // psMesh->setSurfaceColor(glm::vec3(0.255, 0.514, 0.996));
        // psMesh->setEdgeWidth(1.0);
        polyscope::state::userCallback = [&]() { sceneCallback(); };
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
        if (ImGui::Button("advanceOneStep")) 
        {
            simulation.advanceOneStep(static_solve_step++);
            std::vector<glm::vec3> nodes;
            
            for (int i = 0; i < simulation.deformed.rows() / 2; i++)
            {
                TV node_i = simulation.deformed.segment(i * 2, 2);
                nodes.emplace_back(node_i[0], node_i[1], 0);
            }
            
            rod->updateNodePositions(nodes);
            
        }
        if (ImGui::Button("DerivativeTest")) 
        {
            // simulation.testGradient2ndOrderTerm();
            // simulation.testHessian2ndOrderTerm();
            // simulation.testGradientFD();
        }
        if (!animate_modes && run_sim)
        {
            bool finished = simulation.advanceOneStep(static_solve_step++);
            std::vector<glm::vec3> nodes;
            
            for (int i = 0; i < simulation.deformed.rows() / 2; i++)
            {
                TV node_i = simulation.deformed.segment(i * 2, 2);
                nodes.emplace_back(node_i[0], node_i[1], 0);
            }
            
            rod->updateNodePositions(nodes);
            if (finished)
            {
                run_sim = false;
            }
        }
    }
    void run() { polyscope::show(); }

public:
    App(Simulation& sim) : simulation(sim) {}
    ~App() {}
};

#endif
