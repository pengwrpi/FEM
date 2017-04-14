#include "FEA.h"
void FEA::Analysis()
{
    Reorder();
    pMeshIter it = mesh->begin(NSD);
    pMeshEnt e;
    while ((e = mesh->iterate(it)))
    {
        ReadEle(e);
        Ele_Stiff();
        Ele_Load();
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
    SolveEq();
    //clean the "global" variable
    for (unsigned int i = 0; i < NSD; ++i)
    {
        for (unsigned int j = 0; j < NNP; ++j)
            delete[] ID[i][j];
        delete[] ID[i];
    }
    delete[] ID;

    for (unsigned int i = 0; i < NEL; ++i)
    {
        delete[] face_index[i];
        if (num_edge[i] == 0)
            delete[] traction[i][0];
        for (unsigned int j = 0; j < num_edge[i]; ++j)
            delete[] traction[i][j];
        delete[] traction[i];
    }
    delete[] num_edge;
    delete[] face_index;
    delete[] traction;
}
