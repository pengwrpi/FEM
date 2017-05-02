#include "FEA.h"

//This function calculate element load vector contributed by body force

void FEA::Ele_body()
{
    int NINT = numberofintegrationpoints(NEN, NSD);
    double** xilist;
    double* w;
    integrationpoints(NEN, NINT, xilist, NSD);
    integrationweights(NEN, NINT, w, NSD);

    for (unsigned int intpt = 0; intpt < NINT; ++intpt)
    {
        double xi[NSD];
        for (unsigned int i = 0; i < NSD; ++i)
            xi[i] = xilist[i][intpt];
        //compute shape functions and derivatives
        double* N;
        shapefunctions(NEN, xi, N, NSD);
        double** dNdxi;
        shapefunctionderivs(NEN, xi, dNdxi, NSD);
        //compute jacobian matrix & its determinant
        double dxdxi[NSD][NSD];
        for (unsigned int i = 0; i < NSD; ++i)
        {
            for (unsigned int j = 0; j < NSD; ++j)
            {
                dxdxi[i][j] = 0.0;
                for (unsigned int b = 0; b < NEN; ++b)
                    dxdxi[i][j] += coord[i][b] * dNdxi[b][j];
            }
        }
        double dt = dxdxi[0][0] * dxdxi[1][1] - dxdxi[0][1] * dxdxi[1][0]; 
        for (unsigned int a = 0; a < NEN; ++a)
        {
            for (unsigned int i = 0; i < NSD; ++i)
            {
                int p = NSD * a + i;
                FE[p] += N[a] * body_f[i] * w[intpt] * dt;
            }
        }
        for (unsigned int a = 0; a < NEN; ++a)
            delete[] dNdxi[a];
        delete[] dNdxi, N;
    }
    //clean xilist, w
    for (unsigned int i = 0; i < NSD; ++i)
        delete[] xilist[i];
    delete[] xilist;
    delete[] w;
    std::cout << "Ele_body for element " << Ele_ID << " is done" << std::endl;
}
