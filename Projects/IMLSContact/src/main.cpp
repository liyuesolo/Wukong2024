#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/readOBJ.h>
#include "../include/IMLSContact.h"
#include "../include/App.h"

int main()
{
    IMLSContact imls_contact;
    imls_contact.initializeFromFile("../../../Projects/IMLSContact/data/knot.obj");
    // discrete_shell.generateAutodiffRuntime();
    // discrete_shell.setHingeStiffness();

    // for (int j = 80; j < 100; j++)
    // {
    //     for (int d = 0; d < 3; d++)
    //     {
    //         discrete_shell.dirichlet_data[j * 3 + d] = 0;
    //     }
    // }
    
    imls_contact.gravity[1] = 9.8;
    
    App<IMLSContact> app(imls_contact);
    app.initializeScene();
    app.run();
    
    return 0;
}
