#ifndef APP_H
#define APP_H

#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "EoLRodSim.h"

class App
{
public:
    using TV = Vector<T, 3>;
    EoLRodSim& simulation;
    int static_solve_step = 0;

    polyscope::CurveNetwork* rod_network;
    polyscope::PointCloud* rod_vertices;

    T t = 0;

    bool animate_modes = false;
    bool run_sim = false;
    int modes = 0;
    Eigen::MatrixXd eigen_vectors;
    Eigen::VectorXd eigen_values;

    Eigen::MatrixXd meshV;
    Eigen::MatrixXi meshF;

public:
    void initializeScene();

    void sceneCallback();

    void run()
    {
        polyscope::show();
    }

public:
    App(EoLRodSim& sim) : simulation(sim) {}
    ~App() {}
};

#endif