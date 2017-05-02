#include "FEA.h"

//Assemble the global stiffness and load vector
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
                F[P] += FE[p];
                int q = 0;
                for (unsigned int b = 0; b < NEN; ++b)
                {
                    for (unsigned int j = 0; j < NSD; ++j)
                    {
                        if (LM[j][b][0] == 0)//if dof
                        {
                            int Q = LM[j][b][1];
                            //store the information needed for global stiffness matrix
                            Coeff.push_back(T(P, Q, KE[p][q]));
                        }
                        if (LM[j][b][0] == 1) //if dog
                        {
                            int Q = LM[j][b][1];
                            F[P] -= G[Q] * KE[p][q];
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
