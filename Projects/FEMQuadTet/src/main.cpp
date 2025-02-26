#include "../include/FEMQuadTet.h"
#include "../include/App.h"

int main()
{

	FEMQuadTet simulation;
	simulation.initializeFromFile("../../../Projects/FEMQuadTet/data/beam.msh");

	App<FEMQuadTet> app(simulation);
	app.initializeScene();
	app.run();

	return 0;
}