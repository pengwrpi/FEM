#include "FEA.h"
    
int FEA::numberofintegrationpoints(const int& NEN_, const int& NSD_)
{
    if (NSD_ == 2)
    {
        switch (NEN_)
        {
            case 3:
                return 1;
                break;
            case 4:
                return 4;
                break;
            case 6:
                return 3;
                break;
            case 8:
                return 9;
                break;
        }
    } else if (NEN_ == 1) {
        return NEN_;
    }
}

void FEA::integrationpoints(const int& NEN_, const int& NINT, double**& xi, const int& NSD_)
{
    xi = new double*[NSD_];
    for (unsigned int i = 0; i < NSD_; ++i)
        xi[i] = new double[NINT]();
    if (NSD_ == 2)
    {
        if (NEN_ == 3 || NEN_ == 6)
        {
            if (NINT == 1)
            {
                xi[0][0] = 1.0/3.0;
                xi[1][0] = 1.0/3.0;
            }
            if (NINT == 3)
            {
                xi[0][0] = 0.6;
                xi[1][0] = 0.2;
                xi[0][1] = 0.2;
                xi[1][1] = 0.6;
                xi[0][2] = 0.2;
                xi[1][2] = 0.2;
            }
            if (NINT == 4)
            {
                xi[0][0] = 1.0/3.0;
                xi[1][0] = 1.0/3.0;
                xi[0][1] = 0.6;
                xi[1][1] = 0.2;
                xi[0][2] = 0.2;
                xi[1][2] = 0.6;
                xi[0][3] = 0.2;
                xi[1][3] = 0.2;
            }
        }
        if (NEN_ == 4 || NEN_ == 8)
        {
            if (NINT == 1)
            {
                xi[0][0] = 0.0;
                xi[1][0] = 0.0;
            }
            if (NINT == 4)
            {
                xi[0][0] = -0.5773502692; 
                xi[1][0] = xi[0][0];
                xi[0][1] = -xi[0][0];
                xi[1][1] = xi[0][0];
                xi[0][2] = xi[0][0];
                xi[1][2] = -xi[0][0];
                xi[0][3] = -xi[0][0];
                xi[1][3] = -xi[0][0];
            }
            if (NINT == 9)
            {
                xi[0][0] = -0.7745966692;
                xi[1][0] = xi[0][0];
                xi[0][1] = 0.0; 
                xi[1][1] = xi[0][0];
                xi[0][2] = -xi[0][0]; 
                xi[1][2] = xi[0][0];
                xi[0][3] = xi[0][0];
                xi[1][3] = 0.0; 
                xi[0][4] = 0.0;
                xi[1][4] = 0.0;
                xi[0][5] = -xi[0][0];
                xi[1][5] = 0.0; 
                xi[0][6] = xi[0][0];
                xi[1][6] = -xi[0][0]; 
                xi[0][7] = 0.0; 
                xi[1][7] = -xi[0][0];
                xi[0][8] = -xi[0][0]; 
                xi[1][8] = -xi[0][0]; 
            } 
        }   
    } else if (NSD_ == 1) {
        switch (NINT)
        {
            case 1:
                xi[0][0] = 0.0;
                break;
            case 2:
                xi[0][0] = -0.5773502692;
                xi[0][1] = -xi[0][0];
                break;
            case 3:
                xi[0][0] = -0.7745966692;
                xi[0][1] = 0.0;
                xi[0][2] = -xi[0][0];
                break;
        }
    }
}

void FEA::integrationweights(const int& NEN_, const int& NINT, double*& w, const int& NSD_)
{
    w = new double[NINT];
    if (NSD_ == 2)
    {
        if (NEN_ == 3 || NEN_ == 6)
        {
            if (NINT == 1)
            {
                w[0] = 0.5;
            }
            if (NINT == 3)
            {
                w[0] = 1.0/6.0;
                w[1] = 1.0/6.0;
                w[2] = 1.0/6.0;
            }
            if (NINT == 4)
            {
                w[0] = -27.0/96.0;
                w[1] = 25.0/96.0;
                w[2] = 25.0/96.0;
                w[3] = 25.0/96.0;
            }
        }
        if (NEN_ == 4 || NEN_ == 8)
        {
            if (NINT == 1)
            {
                w[0] = 4.0;
            }
            if (NINT == 4)
            {
                w[0] = 1.0;
                w[1] = 1.0;
                w[2] = 1.0;
                w[3] = 1.0;
            }
            if (NINT == 9)
            {
                double w1D[3];
                w1D[0] = 0.555555555;
                w1D[1] = 0.888888888;
                w1D[2] = 0.555555555;
                for (unsigned int j = 0; j < 3; ++j)
                {
                    for (unsigned int i = 0; i < 3; ++i)
                    {
                        int n = 3 * j + i;
                        w[n] = w1D[i] * w1D[j];
                    }
                }
            }
        }
    } else if (NSD_ == 1) {
        switch(NINT)
        {
            case 1:
                w[0] = 2.0;
                break;
            case 2:
                w[0] = 1.0;
                w[1] = 1.0;
                break;
            case 3:
                w[0] = 0.555555555555;
                w[1] = 0.888888888888;
                w[2] = 0.555555555555;
                break;
        }
    }
}

void FEA::shapefunctions(const int NEN_, const double* xi, double*& N, const int& NSD_)
{
    N = new double[NEN_];
    if (NSD_ == 2)
    {
        switch (NEN_)
        {
            case 3:
                N[0] = xi[0];
                N[1] = xi[1];
                N[2] = 1.0 - xi[0] - xi[1];
                break;
            case 6:
                {   
                    double x3 = 1.0 - xi[0] - xi[1];
                    N[0] = (2.0 * xi[0] - 1.0) * xi[0]; 
                    N[1] = (2.0 * xi[1] - 1.0) * xi[1]; 
                    N[2] = (2.0 * x3 - 1.0) * x3;
                    N[3] = 4.0 * xi[0] * xi[1];
                    N[4] = 4.0 * xi[1] * x3;
                    N[5] = 4.0 * xi[0] * x3;
                    break;
                }
            case 4:
                N[0] = 0.25 * (1.0 - xi[0]) * (1.0 - xi[1]);
                N[1] = 0.25 * (1.0 + xi[0]) * (1.0 - xi[1]);
                N[2] = 0.25 * (1.0 + xi[0]) * (1.0 + xi[1]);
                N[3] = 0.25 * (1.0 - xi[0]) * (1.0 + xi[1]);
                break;
            case 8:
                N[0] = -0.25 * (1.0 - xi[0]) * (1.0- xi[1]) * 
                        (1.0 + xi[0] + xi[1]);
                N[1] = 0.25 * (1.0 + xi[0]) * (1.0- xi[1]) * 
                        (xi[0] - xi[1] - 1.0);
                N[2] = 0.25 * (1.0 + xi[0]) * (1.0+ xi[1]) * 
                        (xi[0] + xi[1] - 1.0);
                N[3] = 0.25 * (1.0 - xi[0]) * (1.0+ xi[1]) * 
                        (xi[1] - xi[0] - 1.0);
                N[4] = 0.5 * (1.0 - xi[0] * xi[0]) * (1.0 - xi[1]);
                N[5] = 0.5 * (1.0 + xi[0]) * (1 - xi[1] * xi[1]);
                N[6] = 0.5 * (1.0 - xi[0] * xi[0]) * (1.0 + xi[1]);
                N[7] = 0.5 * (1.0 - xi[0]) * (1 - xi[1] * xi[1]);
                break;
        }
    } else if (NSD_ == 1) {
        switch(NEN_)
        {
            case 2:
                N[0] = 0.5 * (1.0 + xi[0]);
                N[1] = 0.5 * (1.0 - xi[0]);
                break;
            case 3:
                N[0] = -0.5 * xi[0] * (1.0 - xi[0]);
                N[1] = 0.5 * xi[0] * (1.0 + xi[0]);
                N[2] = (1.0 - xi[0]) * (1.0 + xi[0]);
                break;
        }
    }
}

void FEA::shapefunctionderivs(const int& NEN_, const double* xi, double**& dNdxi, const int& NSD_)
{
    dNdxi = new double*[NEN_];
    for (unsigned int i = 0; i < NEN_; ++i)
    {
        dNdxi[i] = new double[NSD_];
    }
    if (NSD_ == 2)
    {
        switch (NEN_)
        {
            case 3:
                dNdxi[0][0] = 1.0;
                dNdxi[0][1] = 0.0;
                dNdxi[1][0] = 0.0;
                dNdxi[1][1] = 1.0;
                dNdxi[2][0] = -1.0;
                dNdxi[2][1] = -1.0;
                break;
            case 6:
                {
                    double x3 = 1 - xi[0] - xi[1];
                    dNdxi[0][0] = 4.0 * xi[0] - 1.0;
                    dNdxi[0][1] = 0.0;
                    dNdxi[1][0] = 0.0;
                    dNdxi[1][1] = 4.0 * xi[1] - 1.0;
                    dNdxi[2][0] = -(4.0 * x3 - 1.0); 
                    dNdxi[2][1] = -(4.0 * x3 - 1.0); 
                    dNdxi[3][0] = 4.0 * xi[1];
                    dNdxi[3][1] = 4.0 * xi[0];
                    dNdxi[4][0] = -4.0 * xi[1]; 
                    dNdxi[4][1] = -4.0 * xi[0];
                    dNdxi[5][0] = 4.0 * x3 - 4.0 * xi[0];
                    dNdxi[5][1] = 4.0 * x3 - 4.0 * xi[1];
                    break;
                }
            case 4:
                dNdxi[0][0] = -0.25 * (1.0 - xi[1]); 
                dNdxi[0][1] = -0.25 * (1.0 - xi[0]);
                dNdxi[1][0] = 0.25 * (1.0 - xi[1]); 
                dNdxi[1][1] = -0.25 * (1.0 + xi[0]);
                dNdxi[2][0] = 0.25 * (1.0 + xi[1]); 
                dNdxi[2][1] = 0.25 * (1.0 + xi[0]);
                dNdxi[3][0] = -0.25 * (1.0 + xi[1]); 
                dNdxi[3][1] = 0.25 * (1.0 - xi[0]);
                break;
            case 8:
                dNdxi[0][0] = 0.25 * (1 - xi[1]) * (2 * xi[0] + xi[1]); 
                dNdxi[0][1] = 0.25 * (1 - xi[0]) * (xi[0] + 2 * xi[1]); 
                dNdxi[1][0] = 0.25 * (1 - xi[1]) * (2 * xi[0] - xi[1]); 
                dNdxi[1][1] = 0.25 * (1 + xi[0]) * (2 * xi[1] - xi[0]); 
                dNdxi[2][0] = 0.25 * (1 + xi[1]) * (2 * xi[0] + xi[1]); 
                dNdxi[2][1] = 0.25 * (1 + xi[0]) * (2 * xi[1] + xi[0]); 
                dNdxi[3][0] = 0.25 * (1 + xi[1]) * (2 * xi[0] - xi[1]); 
                dNdxi[3][1] = 0.25 * (1 - xi[0]) * (2 * xi[1] - xi[0]); 
                dNdxi[4][0] = -xi[0] * (1 - xi[1]);
                dNdxi[4][1] = -0.5 * (1 - xi[0] * xi[0]);
                dNdxi[5][0] = 0.5 * (1 - xi[1] * xi[1]);
                dNdxi[5][1] = -(1 + xi[0]) * xi[1];
                dNdxi[6][0] = -xi[0] * (1 + xi[1]);
                dNdxi[6][1] = 0.5 * (1- xi[0] * xi[0]);
                dNdxi[7][0] = -0.5 * (1-xi[1] * xi[1]);
                dNdxi[7][1] = -(1 - xi[0]) * xi[1];
                break;
        }
    } else {
        switch (NEN_)
        {
            case 2:
                dNdxi[0][0] = 0.5;
                dNdxi[1][0] = -0.5;
                break;
            case 3:
                dNdxi[0][0] = -0.5 + xi[0];
                dNdxi[1][0] = 0.5 + xi[0];
                dNdxi[2][0] = -2.0 * xi[0];
                break;
        }
    }
}

void FEA::materialstiffness(const int& NSD_, const double* materialprops_, double****& C)
{
    C = new double***[NSD_];
    for (unsigned int i = 0; i < NSD_; ++i)
    {
        C[i] = new double**[NSD_];
        for (unsigned int j = 0; j < NSD_; ++j)
        {
            C[i][j] = new double*[NSD_];
            for (unsigned int k = 0; k < NSD_; ++k)
            {
                C[i][j][k] = new double[NSD_]();
            }
        }
    }
    double mu = materialprops_[0];
    double nu = materialprops_[1];
    int planestrain = materialprops_[2];
    for (unsigned int i = 0; i < NSD_; ++i)
    {
        for (unsigned int j = 0; j < NSD_; ++j)
        {
            for (unsigned int k = 0; k < NSD_; ++k)
            {
                for (unsigned int l = 0; l < NSD_; ++l)
                {
                    if (planestrain == 1)
                    {
                        if (i == j && k == l)
                        {
                            C[i][j][k][l] += 2 * mu * nu / (1.0 - 2 * nu);
                        }
                    } else {
                        if (i == j && k == l)
                        {
                            C[i][j][k][l] += 2 * mu * nu / (1 - nu);
                        }
                    }
                    if (i == l && k == j)
                        C[i][j][k][l] += mu;
                    if (i == k && j == l)
                        C[i][j][k][l] += mu;
                }
            }
        }
    }
}

int FEA::nfacenodes(const int& NSD_, const int& NEN_)
{
    if (NSD_ != 2)
        std::cerr << "Only able to handle 2D case now" << std::endl;
    else
    {
        switch(NEN_)
        {
            case 3:
                return 2;
                break;
            case 4:
                return 2;
                break;
            case 6:
                return 3;
                break;
            case 8:
                return 3;
                break;
        }
    }
}

void FEA::facenodes(const int& NSD_, const int& NEN_, 
                    const int& face_index, 
                    int*& node_list/*output*/ )
{
    if (NSD_ != 2)
        std::cerr << "Only able to handle 2D case now" << std::endl;
    else
    {
        int i3[3] = {1, 2, 0};
        int i4[4] = {1, 2, 3, 0};
        node_list = new int[nfacenodes(NSD_, NEN_)];
        switch (NEN_)
        {
            case 3:
                node_list[0] = face_index;
                node_list[1] = i3[face_index];
                break;
            case 6:
                node_list[0] = face_index;
                node_list[1] = i3[face_index];
                node_list[2] = face_index + 3;
                break;
            case 4:
                node_list[0] = face_index;
                node_list[1] = i4[face_index];
                break;
            case 8:
                node_list[0] = face_index;
                node_list[1] = i4[face_index];
                node_list[2] = face_index + 4;
                break;
        }
    }
}
