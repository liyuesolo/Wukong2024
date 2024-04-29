#include "../include/IsohedralTiling.h"
#include "../include/App.h"

int main()
{
    IsohedralTiling tiling;
    tiling.generateOnePerodicUnit();
    App<IsohedralTiling> app(tiling);
    app.initializeScene();
    app.run();
    return 0;
}