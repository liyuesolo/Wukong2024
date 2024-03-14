#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "../include/MMC2D.h"
#include "../include/MMC3D.h"
#include "../include/App.h"

int main()
{
    MMC3D simulation;
    simulation.initialize();
    App<MMC3D> app(simulation);
    app.initializeScene();
    app.run();
    
    return 0;
}