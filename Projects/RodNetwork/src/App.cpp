#include "../include/App.h"

void App::sceneCallback()
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
        simulation.advanceOneStep(static_solve_step++);
        std::vector<glm::vec3> nodes;
        for(auto& rod : simulation.rods)
        {
            for (int idx : rod->indices)
            {
                TV node_i;
                rod->x(idx, node_i);
                nodes.emplace_back(node_i[0], node_i[1], node_i[2]);
            }
        }
        rod_network->updateNodePositions(nodes);
        rod_vertices->updatePointPositions(nodes);
        std::cout << simulation.deformed_states.transpose() << std::endl;
        
    }
    if (ImGui::Button("DerivativeTest")) 
    {
        // simulation.checkTotalGradientScale(true);
    }
    if (!animate_modes && run_sim)
    {
        bool finished = simulation.advanceOneStep(static_solve_step++);
        std::vector<glm::vec3> nodes;
        for(auto& rod : simulation.rods)
        {
            for (int idx : rod->indices)
            {
                TV node_i;
                rod->x(idx, node_i);
                nodes.emplace_back(node_i[0], node_i[1], node_i[2]);
            }
        }
        rod_network->updateNodePositions(nodes);
        rod_vertices->updatePointPositions(nodes);
        if (finished)
        {
            run_sim = false;
        }
    }
}

void App::initializeScene()
{
    polyscope::options::autocenterStructures = true;
    polyscope::view::windowWidth = 3000;
    polyscope::view::windowHeight = 2000;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;

    // Initialize polyscope
    polyscope::init();
    // std::cout << simulation.Rods.size() << std::endl;
    std::vector<glm::vec3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    int global_cnt = 0;
    for(auto& rod : simulation.rods)
    {
        for (int idx : rod->indices)
        {
            TV node_i;
            rod->x(idx, node_i);
            nodes.push_back(glm::vec3(node_i[0], node_i[1], node_i[2]));
        }
        for (int i = 0; i < rod->indices.size()-1; i++)
        {
            edges.push_back({size_t(global_cnt + i), size_t(global_cnt + i + 1)});
        }
        global_cnt += rod->indices.size();
    }
    rod_network = polyscope::registerCurveNetwork("rod network", nodes, edges);
    rod_vertices = polyscope::registerPointCloud("nodes", nodes);
    rod_vertices->setPointRadius(0.007);
    
    rod_network->setColor(glm::vec3(0.255, 0.514, 0.996));
    // psMesh->setEdgeWidth(1.0);

    polyscope::state::userCallback = [&](){ sceneCallback(); };
}