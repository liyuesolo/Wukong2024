#include "../include/DiscreteShellMacro.h"
#include "../include/NeuralMaterialModel.h"
#include "../include/App.h"

int main()
{
    std::string base_folder = "/home/yueli/Documents/ETH/WuKong2024/Projects/";

    cppflow::model tf_model(base_folder + "NeuralMaterialModel/python/Models/10");
    NeuralMaterialModel neural_model(tf_model);
    // neural_model.queryNetworkDerivatives();

    DiscreteShellMacro<NeuralMaterialModel> discrete_shell(neural_model);
    discrete_shell.initializeFromFile(base_folder + "/DiscreteShell/data/grid.obj");
    for (int j = 80; j < 100; j++)
    {
        for (int d = 0; d < 3; d++)
        {
            discrete_shell.dirichlet_data[j * 3 + d] = 0;
        }
    }
    
    discrete_shell.gravity[1] = 0.098;
    
    App<DiscreteShellMacro<NeuralMaterialModel>> app(discrete_shell);
    app.initializeScene();
    app.run();
    return 0;
}