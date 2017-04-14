#include "FEA.h"

void FEA::Ele_Load()
{
    int nfnodes = nfacenodes(NSD, NEN); 
    int NSD_face = NSD - 1;
    int NINT = numberofintegrationpoints(NSD_face, nfnodes);
    double xi[NSD_face];
    double dxdxi[NSD][NSD_face];
    FE = new double[NSD*NEN]();
    double** xilist;
    double* w;
    integrationpoints(nfnodes, NINT, xilist, NSD_face);
    integrationweights(nfnodes, NINT, w, NSD_face);
    //read(num_edge, face_index, traction);
    //traction should be num_edge * NSD
    //num_edge should be how many edges feel the load
    //face_index should be 1D array of face index, with length of num_edge
    for(unsigned int edgeinx = 0; edgeinx < num_edge[Ele_ID]; ++edgeinx)
    {
        int* node_list;
        std::cout << edgeinx << " " << Ele_ID << std::endl;
        facenodes(NSD, NEN, face_index[Ele_ID][edgeinx], node_list); 
        for (unsigned int intpt = 0; intpt < NINT; ++intpt)
        {
            for (unsigned int i = 0; i < NSD_face; ++i)
                xi[i] = xilist[i][intpt];
            double* N;
            shapefunctions(nfnodes, xi, N, NSD_face);
            double** dNdxi;
            shapefunctionderivs(nfnodes, xi, dNdxi, NSD_face);
            //compute jacobian matrix & its determinant
            for (unsigned int i = 0; i < NSD; ++i)
            {
                for (unsigned int j = 0; j < NSD_face; ++j)
                {
                    dxdxi[i][j] = 0.0;
                    for (unsigned int a = 0; a < nfnodes; ++a)
                        dxdxi[i][j] += coord[i][a] * dNdxi[a][j];
                }
            }
            double dt = sqrt(dxdxi[0][0] * dxdxi[0][0] +
                             dxdxi[1][0] * dxdxi[1][0]);
            for (unsigned int a = 0; a < nfnodes; ++a)
            {
                for (unsigned int i = 0; i < NSD; ++i)
                {
                    int row = NSD * node_list[a] + i;
                    FE[row] += N[a] * traction[Ele_ID][edgeinx][i] * w[intpt] * dt;
                }
            }
            for (unsigned int a = 0; a < nfnodes; ++a)
                delete[] dNdxi[a];
            delete[] dNdxi, N;
        }
        delete[] node_list;
    }
    std::cout << "Ele_Load for element " << Ele_ID << " is done" << std::endl;
    //clean xilist, w
    for (unsigned int i = 0; i < NSD_face; ++i)
        delete[] xilist[i];
    delete[] xilist;
    delete[] w;
} 
