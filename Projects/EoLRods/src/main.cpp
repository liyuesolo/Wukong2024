#include "../include/EoLRodSim.h"
#include "../include/Scene.h"
#include "../include/App.h"
int main()
{
    EoLRodSim sim;
    Scene scene(sim);
    scene.buildInterlockingSquareScene(32);
    
    App<EoLRodSim> app(sim);
    app.initializeScene();
    app.run();
    return 0;
}