
#include "../include/NeuralMaterialModel.h"
#include "../include/Timer.h"

void NeuralMaterialModel::queryNetworkDerivatives()
{
    std::vector<double> green_strain_batch = {
        -0.03333, 0.0, 0.03333, -0.03333, 0.0, 0.0, 0.0, 0.0, 0.0,
        -0.03333, 0.0, 0.03333, -0.03333, 0.0, 0.0, 0.0, 0.0, 0.0
    };

    // std::vector<double> green_strain_batch(18 * 1, 1.0);
    
    // cppflow::tensor strain_tf_tensor(green_strain, {1, 3});
    cppflow::tensor strain_tf_tensor(green_strain_batch, {1, 18});

    cppflow::model model("../../../Projects/NeuralMaterialModel/python/Models/10");
    for (int i = 0; i < 10; i++)
    {
        START_TIMING(prediction)
        auto output = model(
            {
                {"serving_default_input_1:0", strain_tf_tensor}
            },
            {
                "StatefulPartitionedCall:0", 
                "StatefulPartitionedCall:1",
                "StatefulPartitionedCall:2"
            }
        );
        FINISH_TIMING_PRINT(prediction)
        std::cout << "energy " << output[0] << std::endl;
        std::cout << "grad " << output[1] << std::endl;
        std::cout << "hessian " << output[2] << std::endl;
        std::getchar();
    }
    
    
    
}