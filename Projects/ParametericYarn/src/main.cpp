#include "../include/ParametericYarn.h"
#include "../include/App.h"

int main()
{

	ParametericYarn parameteric_yarn;
	
	App<ParametericYarn> app(parameteric_yarn);
	app.initializeScene();
	app.run();

	return 0;
}
