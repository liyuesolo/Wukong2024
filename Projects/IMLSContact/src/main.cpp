#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/readOBJ.h>
#include "../include/IMLSContact.h"
#include "../include/App.h"

int main()
{
    IMLSContact simulation;
    simulation.initializeSingleKnot();
    simulation.checkTotalHessianScale(true);
    App<IMLSContact> app(simulation);
    app.initializeScene();
    app.run();
    
    return 0;
}