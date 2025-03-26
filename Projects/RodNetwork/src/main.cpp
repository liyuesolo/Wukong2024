#include "../include/App.h"
#include "../include/RodNetwork.h"

int main()
{
    RodNetwork rod_network;
    rod_network.initializeFromFile(
        "../../../Projects/RodNetwork/data/rods.txt");
    rod_network.verbose = true;
    App app(rod_network);
    app.initializeScene();
    app.run();

    return 0;
}