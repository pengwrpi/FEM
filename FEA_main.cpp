#include <pumi.h>
#include <apf.h>
#include <PCU.h>
#include "FEA.h"

int main(int argc, char** argv)
{
    
    if (argc != 4)
        std::cerr << "wrong command line argument" << std::endl;
    
    MPI_Init(&argc, &argv);
    pumi_start();

    FEA project(argv[1], argv[2], argv[3]); 
    project.Analysis();

//        clean LM, IEN, coord, KE, materialsprop,traction, FE;


  //  clean ID;
    return 0;
}
