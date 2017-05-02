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

//FEA is the class which is responsible to do the FEM analysis

class FEA {
    private:
        //member variables
        std::string ele_type; //linear or quadratic
        pGeom geom;//PUMI geom model
        pMesh mesh;//PUMI mesh
        pNumbering numbers;//PUMI numbering after reordered
        int NNP;  //number of nodes in total
        int NEL;  //number of elements in total
        int NSD;  //number of space dimensions (2D for this code, so 2)
        int*** ID;//store the element information (dof or dog)
        int*** LM;//location matrix
        int NDOF; //number of dof in global matrices
        int NDOG; //number of dog associated with non-zero
        std::vector<double> G; //global prescibed displacement vector
        Eig_V F;//global load vector
        std::vector<T> Coeff;//variable used for assembling the global K matrix
        Eig_V d;//global displacement vector
        int NEN; //number of nodes within a specific element
        int* IEN;//global numbering
        double** coord; //coords of nodes within an element NSD * NEN
        double* materialprops;//materials property parameter
        double** KE;//element stiffness matrix
        double* FE;//element load vector
        int Ele_ID;//element ID
        double body_f[3];//global body force vector
        int num_ele_loaded;//number of element loaded by traction
        int* num_edge;//number of edges loaded by traction within an element
        int** face_index;//index of face loaded by traction within an element
        double*** traction;//traction applied to each element

        //private shapefunction and integration rules
        
        //return number of integration points for different elements
        int numberofintegrationpoints(const int& NEN_, const int& NSD_);
        //return integration points coords
        void integrationpoints(const int& NEN_, const int& NINT_, double**& xi_, const int& NSD_);
        //return integration points weight
        void integrationweights(const int& NEN_, const int& NINT_, double*& w_, const int& NSD_);
        //return shape functions
        void shapefunctions(const int NEN_, const double* xi_, double*& N_, const int& NSD_);
        //return shape function derivatives
        void shapefunctionderivs(const int& NEN_, const double* xi_, double**& dNdxi_, const int& NSD_);
        //return stiffness C[i][j][k][l]
        void materialstiffness(const int& NSD_, const double* materialprops_, double****& C_);
        //return number of nodes on each edge that is loaded by traction
        int nfacenodes(const int& NSD_, const int& NEN_);
        //return the list of nodes on edge that is loaded by traction
        void facenodes(const int& NSD_, const int& NEN_, 
                       const int& face_index_, 
                       int*& node_list_/*output*/ );


        //private memberfunction
        
        //reorder dof to minimize the solving expense
        void Reorder();
        //read Element information from PUMI
        void ReadEle(const pMeshEnt& e);
        //calculate Element stiffness matrix
        void Ele_Stiff();
        //calculate the element load vector contributed by traction
        void Ele_Traction();
        //calcualte the element load vector contributed by body force
        void Ele_body();
        //assemble load vector and stiffness matrix
        void Assemble();
        //solving the equation
        void SolveEq();
        //calculate and output the displacement on node, strain, stress on each integration points
        void Ele_Stress(std::fstream& fs);
    public:
        //constructor
        FEA(const char* mesh_file, const char* filename_in, const char* filename_load_);
        //The "main" function of the FEA class
        void Analysis();
};
