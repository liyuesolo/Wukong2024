#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/readOBJ.h>
#include "../include/DiscreteShell.h"
#include "../include/App.h"

int main()
{
    DiscreteShell discrete_shell;
    discrete_shell.initializeFromFile("../../../Projects/DiscreteShell/data/grid.obj");
    // discrete_shell.setHingeStiffness();

    for (int j = 80; j < 100; j++)
    {
        for (int d = 0; d < 3; d++)
        {
            discrete_shell.dirichlet_data[j * 3 + d] = 0;
        }
    }
    
    discrete_shell.gravity[1] = 0.098;
    
    App<DiscreteShell> app(discrete_shell);
    app.initializeScene();
    app.run();
    
    return 0;
}