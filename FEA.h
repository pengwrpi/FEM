#include <fstream>
#include <vector>
#include <Eigen/Sparse>
#include <PCU.h>
#include <pumi.h>
#include <string>
#include <new>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::VectorXd Eig_V;
typedef Eigen::Triplet<double> T;

class FEA {
    private:
        //member variables
        std::string ele_type;
        pGeom geom;
        pMesh mesh;
        pNumbering numbers;
        int NNP;  //number of nodes in total
        int NEL;  //number of elements in total
        int NSD;  //number of space dimensions (2D for this code, so 2)
        int*** ID;
        int*** LM;
        int NDOF; //number of dof in global matrices
        int NDOG; //number of dog associated with non-zero
        std::vector<double> G; //global prescibed displacement vector
        Eig_V F;
        std::vector<T> Coeff;
        Eig_V d;
        int NEN; //number of nodes within a specific element
        int* IEN;//global numbering
        double** coord; //coords of nodes within an element NSD * NEN
        double* materialprops;//materials property parameter
        double** KE;
        double* FE;
        int Ele_ID;
        double body_f[3];
        int* num_edge;
        int** face_index;
        double*** traction;
        //private shapefunction
        int numberofintegrationpoints(const int& NEN_, const int& NSD_);
        void integrationpoints(const int& NEN_, const int& NINT_, double**& xi_, const int& NSD_);
        void integrationweights(const int& NEN_, const int& NINT_, double*& w_, const int& NSD_);
        void shapefunctions(const int NEN_, const double* xi_, double*& N_, const int& NSD_);
        void shapefunctionderivs(const int& NEN_, const double* xi_, double**& dNdxi_, const int& NSD_);
        void materialstiffness(const int& NSD_, const double* materialprops_, double****& C_);
        int nfacenodes(const int& NSD_, const int& NEN_);
        void facenodes(const int& NSD_, const int& NEN_, 
                       const int& face_index_, 
                       int*& node_list_/*output*/ );
        //private memberfunction
        void Reorder();
        void ReadEle(const pMeshEnt& e);
        void Ele_Stiff();
        void Ele_Load();
        void Ele_body();
        void Assemble();
        void SolveEq();
        void Ele_Stress();
    public:
        FEA(const char* mesh_file, const char* filename_in, const char* filename_load_);
        void Analysis();
};
