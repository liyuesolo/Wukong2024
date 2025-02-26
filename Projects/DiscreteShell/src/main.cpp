#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/readOBJ.h>
#include "../include/DiscreteShell.h"
#include "../include/App.h"

int main()
{
    DiscreteShell discrete_shell;
    discrete_shell.initializeDynamicExampleScene("../../../Projects/DiscreteShell/data/grid.obj");
    discrete_shell.verbose = true;

    
    App<DiscreteShell> app(discrete_shell);
    app.initializeScene();
    app.run();
    
    return 0;
}