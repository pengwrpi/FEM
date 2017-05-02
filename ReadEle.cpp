#include "FEA.h"

void FEA::ReadEle(const pMeshEnt& e)
{
    //set NEN, IEN, coord by PUMI library function
    coord = new double*[NSD];
    NEN = pumi_ment_getTopo(e) + 1;
    std::vector<pMeshEnt> adj;
    pumi_ment_getAdj(e, 0, adj);
    if (ele_type.compare("quadratic") == 0)
    {
        NEN = NEN * 2; 
        std::vector<pMeshEnt> adj_edge;
        pumi_ment_getAdj(e, 1, adj_edge);
        adj.insert(adj.end(), adj_edge.begin(), adj_edge.end());
    } 
    std::cout << "NEN = " << NEN << std::endl;
    //initialize IEN, coord
    for (unsigned int i = 0; i < NSD; ++i)
        coord[i] = new double[NEN];
    IEN = new int[NEN];
    for (unsigned int j = 0; j < NEN; ++j)
    {
        double tmp_coord[NSD];
        pumi_node_getCoord(adj[j], 0, tmp_coord);
        for (unsigned int k = 0; k < NSD; ++k)
            coord[k][j] = tmp_coord[k];
        IEN[j] = pumi_ment_getNumber(adj[j], numbers, 0, 0);
    }
    Ele_ID = pumi_ment_getID(e);
    std::cout << "element " << Ele_ID << " is read" << std::endl;
    for (unsigned int i = 0; i < NEN; ++i)
    {
        std::cout << "coord node " << i << " = " << coord[0][i] << " " << coord[1][i] << std::endl;
        std::cout << "IEN node " << i << " = " << IEN[i] << std::endl;
    }

}
