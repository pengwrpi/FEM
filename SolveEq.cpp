#include "FEA.h"
#include <Eigen/Dense>

void FEA::SolveEq()
{
    SpMat K(NDOF, NDOF);
    std::cout << "NDOF = " << NDOF << std::endl;
    K.setFromTriplets(Coeff.begin(), Coeff.end());
    Eigen::SimplicialCholesky<SpMat> chol(K);
    d = chol.solve(F);
    std::cout << "load vector is" << std::endl;
    for (unsigned int i = 0; i < NDOF; ++i)
        std::cout << F[i] << std::endl;
    std::cout << Eigen::MatrixXd(K) << std::endl;
#if 1 
    for (unsigned int i = 0; i < NDOF; ++i)
        std::cout << d[i] << " ";
    std::cout << std::endl;
#endif
    double** displacement = new double*[NNP];
    for (unsigned int i = 0; i < NNP; ++i)
        displacement[i] = new double[2];
    for (unsigned int A = 0; A < NNP; ++A)
    {
        for (unsigned int i = 0; i < NSD; ++i)
        {
            int P = ID[i][A][1];
            switch(ID[i][A][0])
            {
                case 0:
                    displacement[A][i] = d[P];
                    break;
                case 1:
                    displacement[A][i] = G[P];
                    break;
                case 2:
                    displacement[A][i] = 0.0;
                    break;
            }
        }
    }
    std::cout << "Node Displacement x y" << std::endl;
    for (unsigned int i = 0; i < NNP; ++i)
        std::cout << i << " " << displacement[i][0] << " " << displacement[i][1] << std::endl;
    for (unsigned int i = 0; i < NNP; ++i)
        delete[] displacement[i];
    delete[] displacement;
}
