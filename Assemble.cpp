#include "FEA.h"

void FEA::Assemble()
{
    int p = 0;
    for (unsigned int a = 0; a < NEN; ++a)
    {
        for (unsigned int i = 0; i < NSD; ++i)
        {
            if (LM[i][a][0] == 0) //if dof
            {
                int P = LM[i][a][1];
                std::cout << "after F[P], P = " << P << std::endl;
                F[P] += FE[p];
                std::cout << "after F[P], P = " << P << "F[P] = " << FE[p] << std::endl;
                int q = 0;
                for (unsigned int b = 0; b < NEN; ++b)
                {
                    for (unsigned int j = 0; j < NSD; ++j)
                    {
                        if (LM[j][b][0] == 0)//if dof
                        {
                            int Q = LM[j][b][1];
                            //K[P][Q] += KE[p][q];
                            std::cout << "pushed back " << P << " " << Q << " " << KE[p][q] << std::endl;
                            Coeff.push_back(T(P, Q, KE[p][q]));
                            std::cout << "pushed back " << P << " " << Q << " " << KE[p][q] << std::endl;
                        }
                        if (LM[j][b][0] == 1) //if dog
                        {
                            int Q = LM[j][b][1];
                            F[P] -= G[Q] * KE[p][q];
                            std::cout << "when F[P], P = " << P << std::endl;
                        }
                        q++;
                    }
                }
            }
            p++;
        }
    }
    std::cout << "Assemble for element " << Ele_ID << " is done" << std::endl;
}
