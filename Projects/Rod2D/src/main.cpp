#include "../include/Rod2D.h"
#include "../include/App.h"

int main()
{

	Rod2D rod2d;
	rod2d.initialize();
	App<Rod2D> app(rod2d);
	app.initializeScene();
	app.run();

	return 0;
}