#include "FEA.h"

void FEA::Ele_Stress()
{
    double**** dsde;
    materialstiffness(NSD, materialprops, dsde);
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

    //read in elementary displacement matrix DE
    double DE[NSD][NEN];
    for (unsigned int a = 0; a < NEN; ++a)
    {
        for (unsigned int i = 0; i < NSD; ++i)
        {
            if (LM[i][a][0] == 0)
            {
                int P = LM[i][a][1];
                DE[i][a] = d[P];
            } else if (LM[i][a][0] == 1) {
                int P = LM[i][a][1];
                DE[i][a] = G[P];
            } else {
                DE[i][a] = 0.0;
            }
        }
    }
    std::cout << "Element " << Ele_ID << " displacement is " << std::endl;
    for (unsigned int a = 0; a < NEN; ++a)
    {
        for (unsigned int i = 0; i < NSD; ++i)
            std::cout << DE[i][a] << " "; 
        std::cout << std::endl;
    }
        
    std::cout << "Element " << Ele_ID << "strain stress" << std::endl;
    int npoints = numberofintegrationpoints(NEN, NSD);
    std::cout << "npoints for element " << Ele_ID << " is " << npoints << std::endl;
    double** xilist;
    integrationpoints(NEN, npoints, xilist, NSD);
    std::cout << "integration points for element " << Ele_ID << std::endl;
    for (unsigned int i = 0; i < npoints; ++i)
        std::cout << xilist[0][i] << " " << xilist[1][i] << std::endl;
    for (unsigned int intpt = 0; intpt < npoints; ++intpt)
    {
        double xi[NSD];
        for (unsigned int i = 0; i < NSD; ++i)
            xi[i] = xilist[i][intpt];
        double* N;
        shapefunctions(NEN, xi, N, NSD);
        double** dNdxi;
        shapefunctionderivs(NEN, xi, dNdxi, NSD);
        std::cout << "dNdxi is calculated" << std::endl;
        for (unsigned int rows = 0; rows < NEN; ++rows)
        {
            for (unsigned int cols = 0; cols < NSD; ++cols)
                std::cout << dNdxi[rows][cols] << " ";
            std::cout << std::endl;
        }
        //compute the coords of integration points
        double x[NSD];
        for (unsigned int i = 0; i < NSD; ++i)
        {
            x[i] = 0.0;
            for (unsigned int a = 0; a < NEN; ++a)
                x[i] += coord[i][a] * N[a];
        }
        //compute jacobian matrix and inverse
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
        std::cout << "dxdxi is calculated" << std::endl;
        for (unsigned int rows = 0; rows < NSD; ++rows)
        {
            for (unsigned int cols = 0; cols < NSD; ++cols)
                std::cout << dxdxi[rows][cols] << " ";
            std::cout << std::endl;
        }
        double dt = dxdxi[0][0] * dxdxi[1][1] - dxdxi[0][1] * dxdxi[1][0];
        std::cout << "dt = " << dt << std::endl;
        double dxidx[NSD][NSD];
        dxidx[0][0] = dxdxi[1][1] / dt;
        dxidx[0][1] = -dxdxi[0][1] / dt;
        dxidx[1][0] = -dxdxi[1][0] / dt;
        dxidx[1][1] = dxdxi[0][0] / dt;
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
        std::cout << "dNdx matrix is calculated" << std::endl;
        for (unsigned int rows = 0; rows < NEN; ++rows)
        {
            for (unsigned int cols = 0; cols < NSD; ++cols)
                std::cout << dNdx[rows][cols] << " ";
            std::cout << std::endl;
        }
        //calculate strain
        double strain[NSD][NSD];
        for (unsigned int i = 0; i < NSD; ++i) 
        {
            for (unsigned int j = 0; j < NSD; ++j)
            {
                strain[i][j] = 0.0;
                for (unsigned int a = 0; a < NEN; ++a)
                    strain[i][j] += 0.5 * (DE[i][a] * dNdx[a][j] +
                                           DE[j][a] * dNdx[a][i]);
            }
        }
        //calculate stress 
        double stress[NSD][NSD];
        for (unsigned i = 0; i < NSD; ++i)
        {
            for (unsigned int j = 0; j < NSD; ++j)
            {
                stress[i][j] = 0.0;
                for (unsigned int k = 0; k < NSD; ++k)
                {
                    for (unsigned int l = 0; l < NSD; ++l)
                        stress[i][j] += dsde[i][j][k][l] * strain[k][l];
                }
            }
        }
        std::cout << "intpt " << intpt << " " << x[0] << " " << x[1] << " ";
        std::cout << strain[0][0] << " " << strain[0][1] << " "
                  << strain[1][0] << " " << strain[1][1] << " "   
                  << stress[0][0] << " " << stress[0][1] << " "
                  << stress[1][0] << " " << stress[1][1] << " " << std::endl;
    }
    for (unsigned int i = 0; i < NSD; ++i)
        delete[] xilist[i];
    delete[] xilist;
    for (unsigned int i = 0; i < NSD; ++i)
    {
        for (unsigned int j = 0; j < NSD; ++j)
        {
            for (unsigned int k = 0; k < NSD; ++k)
                delete[] dsde[i][j][k];
            delete[] dsde[i][j];
        }
        delete[] dsde[i];
    }
    delete[] dsde;
}

