#include "../include/DiscreteShellMacro.h"
#include "../include/NeuralMaterialModel.h"
#include "../include/App.h"

int main()
{
    cppflow::model tf_model("");
    NeuralMaterialModel neural_model(tf_model);


    DiscreteShellMacro<NeuralMaterialModel> discrete_shell(neural_model);
    discrete_shell.initializeFromFile("../../../Projects/DiscreteShell/data/grid.obj");

    
    App<DiscreteShellMacro<NeuralMaterialModel>> app(discrete_shell);
    app.initializeScene();
    app.run();
    return 0;
}