#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/readOBJ.h>
#include "../include/FEM3D.h"
#include "../include/App.h"

int main()
{
    FEM3D simulation;
    simulation.initializeFromFile("../../../Projects/FEM3D/data/beam.msh");
    App<FEM3D> app(simulation);
    app.initializeScene();
    app.run();
    
    return 0;
}