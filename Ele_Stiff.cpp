#include "FEA.h"

void FEA::Ele_Stiff()
/********************************************************************************/
{
    double**** dsde;
    materialstiffness(NSD, materialprops, dsde);    
    std::cout << "dsde is calculated" << std::endl;
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

    std::cout << "LM results for element " << Ele_ID << std::endl;
    for (unsigned int a = 0; a < NEN; ++a)
        std::cout << LM[0][a][0] << " " << LM[0][a][1] << " " <<
            LM[1][a][0] << " " << LM[1][a][1] << std::endl;

    int NEE = NEN * NSD;
    KE = new double*[NEE];
    for (unsigned int i = 0; i < NEE; ++i)
    {
        KE[i] = new double[NEE]();
    }

    int npoints = numberofintegrationpoints(NEN, NSD);
    std::cout << "npoints for element " << Ele_ID << " is " << npoints << std::endl;
    double** xilist;
    integrationpoints(NEN, npoints, xilist, NSD);
    std::cout << "integration points for element " << Ele_ID << std::endl;
    for (unsigned int i = 0; i < npoints; ++i)
        std::cout << xilist[0][i] << " " << xilist[1][i] << std::endl;
    double* w;
    integrationweights(NEN, npoints, w, NSD);
    std::cout << "integration weights for element " << Ele_ID << std::endl;
    for (unsigned int i = 0; i < npoints; ++i)
        std::cout << w[i] << " ";
    std::cout << std::endl;
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
        std::cout << "dNdxi is calculated" << std::endl;
        for (unsigned int rows = 0; rows < NEN; ++rows)
        {
            for (unsigned int cols = 0; cols < NSD; ++cols)
                std::cout << dNdxi[rows][cols] << " ";
            std::cout << std::endl;
        }
            
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

        std::cout << "dxidx matrix is calculated" << std::endl;
        for (unsigned int rows = 0; rows < NSD; ++rows)
        {
            for (unsigned int cols = 0; cols < NSD; ++cols)
                std::cout << dxidx[rows][cols] << " ";
            std::cout << std::endl;
        }
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
        std::cout << "dNdx matrix is calculated" << std::endl;
        for (unsigned int rows = 0; rows < NEN; ++rows)
        {
            for (unsigned int cols = 0; cols < NSD; ++cols)
            {
                std::cout << dNdx[rows][cols] << " ";
            }
            std::cout << std::endl;
        }
        //double**** dsde;
        //materialstiffness(NSD, materialprops, dsde);    
        //std::cout << "dsde is calculated" << std::endl;
#if 0
        for (unsigned int a = 0; a < NSD; ++a)
        {
            for (unsigned int b = 0; b < NSD; ++b)
            {
                for (unsigned int c = 0; c < NSD; ++c)
                {
                    for (unsigned int d = 0; d < NSD; ++d)
                    {
                        std::cout << dsde[c][d][b][a] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
#endif

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
                                KE[col][row] += dsde[i][j][k][l] * dNdx[b][l] * dNdx[a][j] * w[intpt] * dt;
                            }
                        }
                    }
                }
            }
        }
        //clean dNdxi
        for (unsigned int i = 0; i < NEN; ++i)
            delete[] dNdxi[i];
        delete[] dNdxi;
    }
    std::cout << "Ele_Stiff for element " << Ele_ID << " is done" << std::endl;
    std::cout << "The matrix is" << std::endl;
    for (unsigned int a = 0; a < NEN*NSD; ++a)
    {
        for (unsigned int b = 0; b < NEN*NSD; ++b)
        {
            //if (KE[a][b] <= 1e-10)
              //  KE[a][b] = 0.0;
            std::cout << KE[a][b] << " ";
        }
        std::cout << std::endl;
    }
    //clean xilist, w, dsde 
    for (unsigned int i = 0; i < NSD; ++i)
        delete[] xilist[i];
    delete[] xilist;
    delete[] w;
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
