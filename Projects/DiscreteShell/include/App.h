#ifndef APP_H
#define APP_H

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

template <class Simulation>
class App
{
public:
    using TV = Vector<T, 3>;

    Simulation& simulation;
    int static_solve_step = 0;

    polyscope::SurfaceMesh* psMesh;

    polyscope::PointCloud* p0, *p1, *p2, *p3;
    T t = 0;

    bool animate_modes = false;
    bool run_sim = false;
    int modes = 0;
    Eigen::MatrixXd eigen_vectors;
    Eigen::VectorXd eigen_values;

    Eigen::MatrixXd meshV;
    Eigen::MatrixXi meshF;

    int current_hinge_idx = 0;

public:
    void initializeScene()
    {
        polyscope::options::programName = "WuKongSim";
        // polyscope::options::autocenterStructures = true;
        polyscope::view::windowWidth = 3000;
        polyscope::view::windowHeight = 2000;
        polyscope::options::groundPlaneMode =
            polyscope::GroundPlaneMode::None;
        // polyscope::options::groundPlaneHeightFactor = 0.6;
        // polyscope::options::shadowDarkness = 0.4;

        // Initialize polyscope
        polyscope::init();
        vectorToIGLMatrix<T, 3>(simulation.deformed, meshV);
        vectorToIGLMatrix<int, 3>(simulation.faces, meshF);
        psMesh = polyscope::registerSurfaceMesh("surface mesh", meshV, meshF);
        psMesh->setSmoothShade(false);
        psMesh->setSurfaceColor(glm::vec3(0.255, 0.514, 0.996));
        psMesh->setEdgeWidth(1.0);

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
        if (ImGui::Button("Animate Modes"))
        {
            animate_modes = true;
        }
        if (ImGui::Button("Next Modes"))
        {
            modes = (modes + 1) % eigen_values.rows();
        }
        if (ImGui::Button("Compute Linear Modes"))
        {
            simulation.computeLinearModes(eigen_vectors, eigen_values);
        }
        if (ImGui::Button("advanceOneStep"))
        {
            simulation.advanceOneStep(static_solve_step++);
            vectorToIGLMatrix<T, 3>(simulation.deformed, meshV);
            psMesh->updateVertexPositions(meshV);
            // Eigen::MatrixXd meshV;
            // Eigen::MatrixXi meshF;
            // vectorToIGLMatrix<int, 3>(simulation.faces, meshF);
            // polyscope::registerSurfaceMesh("surface mesh", meshV, meshF);
        }
        if (ImGui::Button("check hinge"))
        {
            auto hinge_idx = simulation.hinges.row(current_hinge_idx);
            auto deformed_vertices = simulation.getHingeVtxDeformed(hinge_idx);
            TV x0 = deformed_vertices.row(0);
            TV x1 = deformed_vertices.row(1);
            TV x2 = deformed_vertices.row(2);
            TV x3 = deformed_vertices.row(3);
            std::vector<glm::vec3> nodes_0 = {glm::vec3(x0[0], x0[1], x0[2])};
            p0 = polyscope::registerPointCloud("x0", nodes_0);
            std::vector<glm::vec3> nodes_1 = {glm::vec3(x1[0], x1[1], x1[2])};
            p1 = polyscope::registerPointCloud("x1", nodes_1);
            std::vector<glm::vec3> nodes_2 = {glm::vec3(x2[0], x2[1], x2[2])};
            p2 = polyscope::registerPointCloud("x2", nodes_2);
            std::vector<glm::vec3> nodes_3 = {glm::vec3(x3[0], x3[1], x3[2])};
            p3 = polyscope::registerPointCloud("x3", nodes_3);
        }
        ImGui::SameLine();
        if (ImGui::Button("next hinge"))
        {
            auto hinge_idx = simulation.hinges.row(++current_hinge_idx);
            auto deformed_vertices = simulation.getHingeVtxDeformed(hinge_idx);
            TV x0 = deformed_vertices.row(0);
            TV x1 = deformed_vertices.row(1);
            TV x2 = deformed_vertices.row(2);
            TV x3 = deformed_vertices.row(3);
            std::vector<glm::vec3> nodes_0 = {glm::vec3(x0[0], x0[1], x0[2])};
            p0->updatePointPositions(nodes_0);
            std::vector<glm::vec3> nodes_1 = {glm::vec3(x1[0], x1[1], x1[2])};
            p1->updatePointPositions(nodes_1);
            std::vector<glm::vec3> nodes_2 = {glm::vec3(x2[0], x2[1], x2[2])};
            p2->updatePointPositions(nodes_2);
            std::vector<glm::vec3> nodes_3 = {glm::vec3(x3[0], x3[1], x3[2])};
            p3->updatePointPositions(nodes_3);
        }
        if (ImGui::Button("check derivative"))
        {
            simulation.checkTotalGradient(true);
            simulation.checkTotalHessian(true);
        }
        ImGui::SameLine();
        if (ImGui::Button("check derivative scale"))
        {
            simulation.checkTotalGradientScale(true);
            simulation.checkTotalHessianScale(true);
        }
        if (ImGui::Checkbox("ConsistentMass",
                            &simulation.use_consistent_mass_matrix))
        {
            if (!simulation.use_consistent_mass_matrix)
                simulation.computeMassMatrix();
        }
        if (ImGui::Checkbox("Dynamics", &simulation.dynamics))
        {
            if (simulation.dynamics)
                simulation.initializeDynamicStates();
        }
        if (ImGui::Checkbox("Gravity", &simulation.add_gravity))
        {
        }
        if (ImGui::Checkbox("Verbose", &simulation.verbose))
        {
        }
        if (animate_modes && !run_sim)
        {
            t += 0.1;
            simulation.deformed = simulation.undeformed + simulation.u +
                                  eigen_vectors.col(modes) * std::sin(t);

            Eigen::MatrixXd meshV;
            vectorToIGLMatrix<T, 3>(simulation.deformed, meshV);
            psMesh->updateVertexPositions(meshV);
        }
        if (!animate_modes && run_sim)
        {
            bool finished = simulation.advanceOneStep(static_solve_step++);
            Eigen::MatrixXd meshV;
            vectorToIGLMatrix<T, 3>(simulation.deformed, meshV);
            psMesh->updateVertexPositions(meshV);
            if (finished)
                run_sim = false;
        }
    }

    void run() { polyscope::show(); }

public:
    App(Simulation& sim) : simulation(sim) {}
    ~App() {}
};

#endif