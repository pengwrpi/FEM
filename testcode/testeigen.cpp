#include <iostream>
#include <Eigen/Sparse>
#include <vector>
int main()
{
    std::vector<double> v1 = {2.3, 2.0, 3.0};
    Eigen::VectorXd v2(v1.size());
    for (unsigned int i = 0; i < v1.size(); ++i)
        v2[i] = v1[i];
    std::cout << "v2 = " << v2 << std::endl;



    Eigen::SparseMatrix<double> A(3,3);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletlist;
    for(unsigned int i = 0; i < 3; ++i)
        tripletlist.push_back(T(i,i,v1[i]));

    A.setFromTriplets(tripletlist.begin(), tripletlist.end());
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(A);
    Eigen::VectorXd x = chol.solve(v2);

    std::cout << "x = " << x << std::endl;

    return 0;

}
