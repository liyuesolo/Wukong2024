#include "../include/RodNetwork.h"
#include "../include/App.h"

int main()
{
	RodNetwork rod_network;
	rod_network.initializeFromFile("../../../Projects/RodNetwork/data/rods.txt");
	App app(rod_network);
	app.initializeScene();
	app.run();

	return 0;
}