#include "FEA.h"

void FEA::Setup(int& NNP,  //number of nodes in total
                int& NEL,  //number of elements in total
                int& NSD,  //number of space dimensions (2D for this code, so 2)
                //std::vector<vector<int> >& ID, //destination matrix
                int*** ID, 
                int& NDOF, //number of dof in global matrices
                int& NDOG, //number of dog associated with non-zero
                std::vector<std::vector<double> >& K,//global stiffness matrix
                std::vector<double>& F, //global force vector
                std::vector<double>& G, //global prescibed displacement vector
                std::fstream& read_in)
/********************************************************************************/
{
    NDOF = 0;
    NDOG = 0;
    ID = new int**[NSD];
    for (unsigned int i = 0; i < NSD; ++i)
    {
        ID[i] = new int*[NNP];
        for (unsigned int j = 0; j < NNP; ++j)
        {
            ID[i][j] = new int[2];
        }
    }
    for (unsigned int A = 0; A < NNP; ++A)
    {
        for (unsigned int i = 0; i < NSD; ++i)
        {
            unsigned int temp_type;
            double FG;
            read_in >> temp_type;
            read_in >> FG;
            switch (temp_type)
            {
                case 0://dof
                    ++NDOF;
                    ID[i][A][0] = 0;
                    ID[i][A][1] = NDOF;
                    F.push_back(FG);
                    break;
                case 1://non-zero dog
                    ++NDOG;
                    ID[i][A][0] = 1;
                    ID[i][A][1] = NDOG;
                    G.push_back(FG);
                    break;
                case 2://zero BC 
                    ID[i][A][0] = 2;
                    ID[i][A][1] = 0;
                    break;
            }
        }
    }
    for (unsigned int A = 0; A < NDOF; ++A)
    {
        std::vector<double> temp;
        for (unsigned int B = 0; B < NDOF; ++B)
        {
            temp.push_back(0.0);
        }
        K.push_back(temp);
        temp.clear();
    }
}

void FEA::Ele_Stiff(const int NSD, 
          //const std::vector<vector<std::pair<int, int> > >& ID, 
          const int*** ID,
          //std::vector<vector<std::pair<int, int> > >& LM, 
          int*** LM, 
          int& NEN,
          //std::vector<vector<double> >& KE,
          double** KE)
/********************************************************************************/
{
    read(NEN, IEN, coord, materialprops);
    LM = new int**[NSD];
    for (unsigned int i = 0; i < NSD; ++i)
    {
        LM[i] = new int*[NEN];
        for (unsigned int j = 0; j < NEN; ++j)
        {
            LM[i][j] = new int[2];
        }
    }

    for (unsigned int a = 0; a < NEN; ++a)
    {
        unsigned int A = IEN[a];
        for (unsigned int i = 0; i < NSD; ++i)
        {
            LM[i][a][0] = ID[i][A][0];
            LM[i][a][1] = ID[i][A][1];
        }
    }
    int NEE = NEN * NSD;
    KE = new double*[NEE];
    for (unsigned int i = 0; i < NEE; ++i)
    {
        KE[i] = new double[NEE];
    }

    int npoints = numberofintegrationpoints(NEN, NSD);
    double** xilist;
    integrationpoints(NEN, npoints, xilist, NSD);
    double* w;
    integrationweights(NEN, npoints, w, NSD);
    for (unsigned int intpt = 0; intpt < npoints; ++intpt)
    {
        //compute shape function & derivatives regards to local coord
        double xi[NSD];
        for (unsigned int i = 0; i < NSD; ++i)
        {
            xi[i] = xilist[i][intpt];
        }
        double** dNdxi;
        shapefunctionderivs(NEN, xi, dNdxi, NSD);
        //compute shape jacobian matrix & its determinant
        double dxdxi[NSD][NSD];
        for (unsigned int i = 0; i < NSD; ++i)
        {
            for (unsigned int j = 0; j < NSD; ++j)
            {
                dxdxi[i][j] = 0.0;
                for (unsigned int a = 0; a < NEN; ++a)
                {
                    dxdxi[i][j]+= coord[i][a] * dNdxi[a][j];
                }
            }
        }
        double det_original = dxdxi[0][0] * dxdxi[1][1] - dxdxi[0][1] * dxdxi[1][0];
        double dxidx[NSD][NSD];
        dxidx[0][0] = dxdxi[1][1] / det_original;
        dxidx[0][1] = -dxdxi[0][1] / det_original;
        dxidx[1][0] = -dxdxi[1][0] / det_original;
        dxidx[1][1] = dxdxi[0][0] / det_original;
        double dt = dxidx[0][0] * dxidx[1][1] - dxidx[0][1] * dxidx[1][0];
        //convert shape function derivatives to global coords
        double dNdx[NEN][NSD];
        for (unsigned int a = 0; a < NEN; ++a)
        {
            for (unsigned int i = 0; i < NSD; ++i)
            {
                dNdx[a][i] = 0.0;
                for (unsigned int j = 0; j < NSD; ++j)
                {
                    dNdx[a][i] += dNdxi[a][j] * dxidx[j][i];
                }
            }
        }
        double**** dsde;
        materialstiffness(NSD, materialprops, dsde);    
        for (unsigned int a = 0; a < NEN; ++a)
        {
            for (unsigned int i = 0; i < NSD; ++i)
            {
                for (unsigned int b = 0; b < NEN; ++b)
                {
                    for (unsigned int k = 0; k < NSD; ++k)
                    {
                        int row = NSD * a + i;
                        int col = NSD * b + k;
                        for (unsigned int j = 0; j < NSD; ++j)
                        {
                            for (unsigned int l = 0; l < NSD; ++l)
                            {
                                KE[col][row] += dsde[i][j][k][l] * dNdx[b][l] *
                                    dNdx[a][j] * w[intpt] * dt;
                            }
                        }
                    }
                }
            }
        }
    }
    clean xilist, w, dNdxi, dsde 
}

void Ele_Load(const int& NSD, const int& NEN, 
              int* node_list, double* FE);
{
    int nfnodes = nfacenodes(NSD, NEN); 
    int NINT = numberofintegrationpoints(NSD - 1, nfnodes);
    double xi[NSD - 1];
    double dxdxi[NSD][NSD - 1];
    FE = new double[NSD*NEN];
    double*** xilist;
    double* w;
    integrationpoints(nfnodes, NINT, xilist, NSD - 1);
    integrationweights(nfnodes, NINT, w, NSD - 1);
    read(coord,num_edge, face_index, traction);
    //traction should be num_edge * NSD
    //num_edge should be how many edges feel the load
    //face_index should be 1D array of face index, with length of num_edge
    for(unsigned int edgeinx = 0; edgeinx < num_edge; ++edgeinx)
    {
        int* node_list;
        facenodes(NSD, NEN, face_index[edgeinx], node_list); 
        for (unsigned int intpt = 0; intpt < NINT; ++intpt)
        {
            for (unsigned int i = 0; i < NSD - 1; ++i)
                xi[i] = xilist[i][intpt];
            double* N;
            shapefunctions(nfnodes, xi, N, NSD - 1);
            double** dNdxi;
            shapefunctionderivs(nfnodes, xi, dNdxi, NSD - 1);
            //compute jacobian matrix & its determinant
            for (unsigned int i = 0; i < NSD; ++i)
            {
                for (unsigned int j = 0; j < NSD - 1; ++j)
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
                    FE[row] += N[a] * traction[edgeinx][i] * w[intpt] * dt;
                }
            }
        }
        delete[] node_list, N;
        for (unsigned int a = 0; a < nfnodes; ++a)
            delete[] dNdxi[a];
        delete[] dNdxi;
    }
    clean xilist, w, N, dNdxi, node_list 
} 

void FEA::Assemble(const int& NSD,
                   const int*** LM, 
                   const int& NEN,
                   const double** KE,
                   const double* FE,
                   const std::vector<double>& G,
                   std::vector<std::vector<double> >& K,
                   std::vector<double>& F)
{
    int p = 0;
    for (unsigned int a = 0; a < NEN; ++a)
    {
        for (unsigned int i = 0; i < NSD; ++i)
        {
            p++;
            if (LM[i][a][0] == 0) //if dof
            {
                int P = LM[i][a][1];
                F[P] += FE[p];
                int q = 0;
                for (unsigned int b = 0; b < NEN; ++b)
                {
                    for (unsigned int j = 0; j < NSD; ++j)
                    {
                        q++;
                        if (LM[j][b][0] == 0)//if dof
                        {
                            int Q = LM[j][b][1];
                            K[P][Q] += KE[p][q];
                        }
                        if (LM[j][b][1] == 1) //if dog
                        {
                            int Q = LM[j][b][1];
                            F[P] -= G[Q] * KE[p][q];
                        }
                    }
                }
            }
        }
    }
}


