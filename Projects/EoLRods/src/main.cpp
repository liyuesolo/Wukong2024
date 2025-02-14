#include "../include/EoLRodSim.h"
#include "../include/Scene.h"
#include "../include/App.h"
int main()
{
    EoLRodSim sim;
    Scene scene(sim);
    // scene.buildInterlockingSquareScene(8);
    // scene.buildFullScaleSquareScene(8);
    scene.buildOneCrossScene(32);
    
    App app(sim);
    app.initializeScene();
    app.run();
    return 0;
}