#include "FEA.h"

//This function compute the element load vector contributed by traction
void FEA::Ele_Traction()
{
    std::cout << "Ele_Load for element " << Ele_ID << " started" << std::endl;
    int nfnodes = nfacenodes(NSD, NEN); 
    std::cout << "nfnodes is " << nfnodes << std::endl;
    std::cout << "NEN is " << NEN << std::endl;
    int NSD_face = NSD - 1;
    int NINT = numberofintegrationpoints(nfnodes, NSD_face);
    std::cout << "NINT is " << NINT << std::endl;
    double xi[NSD_face];
    double dxdxi[NSD][NSD_face];
    FE = new double[NSD*NEN]();
    double** xilist;
    double* w;
    integrationpoints(nfnodes, NINT, xilist, NSD_face);
    integrationweights(nfnodes, NINT, w, NSD_face);
    //read num_edge, face_index, traction;
    //traction should be num_edge * NSD
    //num_edge should be how many edges feel the load
    //face_index should be 1D array of face index, with length of num_edge
    for(unsigned int edgeinx = 0; edgeinx < num_edge[Ele_ID]; ++edgeinx)
    {//loop through all the loaded edges
        int* node_list;
        std::cout << "edgeinx for " << Ele_ID << " is " << edgeinx << " now" << std::endl;
        facenodes(NSD, NEN, face_index[Ele_ID][edgeinx], node_list); 
        std::cout << "node_lists are" << std::endl;
        for (unsigned int i = 0; i < nfnodes; ++i)
            std::cout << node_list[i] << " ";
        std::cout << std::endl;
        std::cout << "coords are" << std::endl;
        for (unsigned int i = 0; i < nfnodes; ++i)
            std::cout << coord[0][node_list[i]] << " " << coord[1][node_list[i]] << std::endl;;
        for (unsigned int intpt = 0; intpt < NINT; ++intpt)
        {//loop through all the integration points
            for (unsigned int i = 0; i < NSD_face; ++i)
                xi[i] = xilist[i][intpt];
            std::cout << "xi[0] = " << xi[0] << std::endl;
            //compute shape functions
            double* N;
            shapefunctions(nfnodes, xi, N, NSD_face);
            std::cout << "N are " << std::endl;
            for (unsigned int a = 0; a < nfnodes; ++a)
                std::cout << N[a] << " ";
            std::cout << std::endl;
            //compute shape functions derivatives
            double** dNdxi;
            shapefunctionderivs(nfnodes, xi, dNdxi, NSD_face);
            std::cout << "dNdxi are " << std::endl;
            for (unsigned int a = 0; a < nfnodes; ++a)
                std::cout << dNdxi[a][0] << " ";
            std::cout << std::endl;
            //compute jacobian matrix & its determinant
            for (unsigned int i = 0; i < NSD; ++i)
            {
                for (unsigned int j = 0; j < NSD_face; ++j)
                {
                    dxdxi[i][j] = 0.0;
                    for (unsigned int a = 0; a < nfnodes; ++a)
                        dxdxi[i][j] += coord[i][node_list[a]] * dNdxi[a][j];
                }
            }
            std::cout << "dxdxi are " << std::endl;
            for (unsigned int i = 0; i < NSD ; ++i)
                std::cout << dxdxi[i][0] << " ";
            std::cout << std::endl;

            //compute determinent
            double dt = sqrt(dxdxi[0][0] * dxdxi[0][0] +
                             dxdxi[1][0] * dxdxi[1][0]);
            std::cout << "dt is " << dt << std::endl;
            for (unsigned int a = 0; a < nfnodes; ++a)
            {
                for (unsigned int i = 0; i < NSD; ++i)
                {
                    int row = NSD * node_list[a] + i;
                    FE[row] += N[a] * traction[Ele_ID][edgeinx][i] * w[intpt] * dt;
                    std::cout << "FE[" << row << "]" << " = " << FE[row] << std::endl;
                }
            }
            for (unsigned int a = 0; a < nfnodes; ++a)
                delete[] dNdxi[a];
            delete[] dNdxi, N;
        }
        delete[] node_list;
    }
    std::cout << "Ele_Load for element " << Ele_ID << " is done" << std::endl;
    for (unsigned int i = 0; i < NSD*NEN; ++i)
        std::cout << FE[i] << " ";
    std::cout << std::endl;
    //clean xilist, w
    for (unsigned int i = 0; i < NSD_face; ++i)
        delete[] xilist[i];
    delete[] xilist;
    delete[] w;
} 
