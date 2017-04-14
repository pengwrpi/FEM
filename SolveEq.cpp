#include "FEA.h"

void FEA::SolveEq()
{
    SpMat K(NDOF, NDOF);
    K.setFromTriplets(Coeff.begin(), Coeff.end());
    Eigen::SimplicialCholesky<SpMat> chol(K);
    d = chol.solve(F);
    for (unsigned int i = 0; i < NDOF; ++i)
        std::cout << d[i] << " ";
    std::cout << std::endl;
}
