#include "FEA.h"

//This is the "main" function which actually runs the FEA
void FEA::Analysis()
{
    //firstly reorder
    Reorder();
    pMeshIter it = mesh->begin(NSD);
    pMeshEnt e;
    while ((e = mesh->iterate(it)))
    {//loop through all the elements
     // and calculate element stiffness matrix
     // load vector and assemble it
        ReadEle(e);
        Ele_Stiff();
        Ele_Traction();
        Ele_body();
        Assemble();
        //clean the variable associated with elements
        for (unsigned int i = 0; i < NSD; ++i)
            delete[] coord[i];
        delete[] coord;
        delete[] IEN;
        for (unsigned int i = 0; i < NSD; ++i)
        {
            for (unsigned int j = 0; j < NEN; ++j)
                delete[] LM[i][j];
            delete[] LM[i];
        }
        delete[] LM;
        for (unsigned int i = 0; i < NEN * NSD; ++i)
            delete[] KE[i];
        delete[] KE;
        delete[] FE;
    }
    mesh->end(it);

    //solving the equation
    SolveEq();

    //create the file which stores the displacement 
    //stress and strain at the integration points
    std::fstream resultfs;
    resultfs.open("FEA_result.txt", std::fstream::out);
    it = mesh->begin(NSD);
    while ((e = mesh->iterate(it)))
    {//loop through all the element
        ReadEle(e);
        Ele_Stress(resultfs);
        //clean the variable associated with elements
        for (unsigned int i = 0; i < NSD; ++i)
            delete[] coord[i];
        delete[] coord;
        delete[] IEN;
        for (unsigned int i = 0; i < NSD; ++i)
        {
            for (unsigned int j = 0; j < NEN; ++j)
                delete[] LM[i][j];
            delete[] LM[i];
        }
        delete[] LM;
    }
    mesh->end(it);
    resultfs.close();
    //clean the "global" variable
    for (unsigned int i = 0; i < NSD; ++i)
    {
        for (unsigned int j = 0; j < NNP; ++j)
            delete[] ID[i][j];
        delete[] ID[i];
    }
    delete[] ID;
    if (num_ele_loaded != 0)
    {
        for (unsigned int i = 0; i < NEL; ++i)
        {
            delete[] face_index[i];
            if (num_edge[i] == 0)
                delete[] traction[i][0];
            for (unsigned int j = 0; j < num_edge[i]; ++j)
                delete[] traction[i][j];
            delete[] traction[i];
        }
    }
    delete[] num_edge;
    delete[] face_index;
    delete[] traction;
    delete[] materialprops;
}
