#include "../include/App.h"
#include "../include/Homogenization.h"
#include "../include/DiscreteShellHomogenization.h"

int main()
{
    DiscreteShellHomogenization shell_homogenizor;
    shell_homogenizor.initializeFromFile("../../../Projects/DiscreteShell/data/grid.obj");

    Homogenization<DiscreteShellHomogenization> macro_model(shell_homogenizor);

    App<Homogenization<DiscreteShellHomogenization>> app(macro_model);
    app.visualization_unit = Triangle;
    app.initializeScene();
    app.run();
    return 0.0;
}