#include "../include/App.h"

template<class Simulation>
void App<Simulation>::initializeScene()
{
    polyscope::options::autocenterStructures = true;
    polyscope::view::windowWidth = 1024;
    polyscope::view::windowHeight = 1024;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::options::groundPlaneHeightFactor = 0.; 
    polyscope::options::shadowDarkness = 0.4;
    polyscope::init();

    // Initialize polyscope
    vectorToIGLMatrix<T, 3>(simulation.native_scale_model.deformed, meshV);
    vectorToIGLMatrix<int, 3>(simulation.native_scale_model.faces, meshF);
    psMesh = polyscope::registerSurfaceMesh("surface mesh", meshV, meshF);
    psMesh->setSmoothShade(false);
    psMesh->setSurfaceColor(glm::vec3(0.255, 0.514, 0.996));
    psMesh->setEdgeWidth(1.0);

    polyscope::state::userCallback = [&](){ sceneCallback(); };
}

template class App<Homogenization<DiscreteShellHomogenization>>;