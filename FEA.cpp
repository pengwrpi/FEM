#include "FEA.h"
#include <fstream>
#include <string>
#include <set>

FEA::FEA(const char* mesh_file, const char* filename_in, const char* filename_load_)
                
/********************************************************************************/
{
    geom = pumi_geom_load(NULL, "null");
    mesh = pumi_mesh_load(geom, mesh_file, 1);
    std::cout << "mesh dimension is " << pumi_mesh_getDim(mesh) << std::endl;
    std::fstream read_in;
    read_in.open(filename_in, std::fstream::in);
    read_in >> ele_type;
    if (ele_type.compare("quadratic") == 0)
    {
        pumi_mesh_setShape(mesh, pumi_shape_getSerendipity());
        NNP = pumi_mesh_getNumEnt(mesh, 0) + pumi_mesh_getNumEnt(mesh, 1);
    } else if (ele_type.compare("linear") == 0) {
        NNP = pumi_mesh_getNumEnt(mesh, 0);
    } else {
        std::cerr << "wrong type of element" << std::endl;
    }
    materialprops = new double[3];
    read_in >> materialprops[0];
    read_in >> materialprops[1];
    read_in >> materialprops[2];

    std::cout << "materialsprops read" << std::endl;
    std::cout << materialprops[0] << " " << materialprops[1] << " " << std::endl; 

    int num_nodes_fix;

    read_in >> num_nodes_fix;
    std::set<int> fixed_nodes;
    int dim[num_nodes_fix];
    int types[num_nodes_fix];
    double condition[num_nodes_fix];

    for(unsigned int i = 0; i < num_nodes_fix; ++i)
    {
        int index_temp;
        read_in >> index_temp;
        fixed_nodes.insert(index_temp);
        read_in >> dim[i];
        read_in >> types[i]; 
        read_in >> condition[i]; 
    }

    read_in.close();
    NSD = pumi_mesh_getDim(mesh);
    NEL = pumi_mesh_getNumEnt(mesh, NSD);
    NDOF = 0;
    NDOG = 0;
    ID = new int**[NSD];
    for (unsigned int i = 0; i < NSD; ++i)
    {
        ID[i] = new int*[NNP];
        for (unsigned int j = 0; j < NNP; ++j)
        {
            ID[i][j] = new int[2];
        }
    }
    std::vector<double> F_std;
    int count_fixed = 0;
    for (unsigned int A = 0; A < NNP; ++A)
    {
        for (unsigned int i = 0; i < NSD; ++i)
        {
            int temp_type;
            double FG;
            if (fixed_nodes.find(A) != fixed_nodes.end() &&
                    dim[count_fixed] == i)
            {
                temp_type = types[count_fixed];
                FG = condition[count_fixed];
                count_fixed++;
            } else {
                temp_type = 0;
                FG = 0.0;
            }
            switch (temp_type)
            {
                case 0://dof
                    ID[i][A][0] = 0;
                    ID[i][A][1] = NDOF;
                    ++NDOF;
                    F_std.push_back(FG);
                    break;
                case 1://non-zero dog
                    ID[i][A][0] = 1;
                    ID[i][A][1] = NDOG;
                    ++NDOG;
                    G.push_back(FG);
                    break;
                case 2://zero BC 
                    ID[i][A][0] = 2;
                    ID[i][A][1] = 0;
                    break;
            }
        }
    }
    std::cout << "element info read" << std::endl;
    for (unsigned int i = 0; i < NNP; ++i)
        std::cout << ID[0][i][0] << " " << ID[0][i][1] << " "
            << ID[1][i][0] << " " << ID[1][i][1] << std::endl;
    F = Eigen::VectorXd::Zero(F_std.size());
    std::cout << "Force size " << F_std.size() << std::endl;
    for (unsigned int i = 0; i < F_std.size(); i++)
        F[i] = F_std[i];
    //read in load information
    num_edge = new int[NEL]();
    face_index = new int*[NEL];
    traction = new double**[NEL];
    std::fstream load_in;
    load_in.open(filename_load_, std::fstream::in);

    for (unsigned int i = 0; i < NSD; ++i)
        load_in >> body_f[i];

    std::set<int> ele_loaded;
    load_in >> num_ele_loaded;
    for (unsigned int i = 0; i < num_ele_loaded; ++i)
    {
        int temp_index;
        load_in >> temp_index;
        ele_loaded.insert(temp_index);
    }

    if (!ele_loaded.empty())
    {
        for (unsigned int i = 0; i < NEL; ++i)
        {
            if (ele_loaded.find(i) == ele_loaded.end())
            {
                num_edge[i] = 0;
                face_index[i] = new int[1];
                traction[i] = new double*[1];
                traction[i][0] = new double[NSD];
                continue;
            }
            int count;
            load_in >> count;
            num_edge[i] = count;
            face_index[i] = new int[count];
            traction[i] = new double*[count];
            for (unsigned int j = 0; j < count; ++j)
                load_in >> face_index[i][j];
            for (unsigned int j = 0; j < count; ++j)
            {
                traction[i][j] = new double[NSD];
                for (unsigned int k = 0; k < NSD; ++k)
                    load_in >> traction[i][j][k];
            }
        }
    }
    load_in.close();
    for (unsigned int i = 0; i < NEL; ++i)
    {
        if (num_edge[i] != 0)
        {
            std::cout << num_edge[i] << " ";
            for (unsigned int j = 0; j < num_edge[i]; ++j)
            {
                std::cout << face_index[i][j] << " ";
                std::cout << traction[i][j][0] << " " << traction[i][j][1] << std::endl;
            }
        }
    }
    std::cout << "load info read" << std::endl;
    
}
