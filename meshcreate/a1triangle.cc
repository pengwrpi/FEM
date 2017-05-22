#include <PCU.h>
#include <pumi.h>
#include <iostream>
#include "createmesh.h"

int main(int argc, char** argv)
{
  if (argc != 3)
      std::cerr << "wrong argument" << std::endl;
  MPI_Init(&argc,&argv);
  pumi_start();
  pGeom geom = pumi_geom_load(NULL, "null");
  pMesh mesh = pumi_mesh_create(geom, 2);

  int x = 1;
  int y = 1;
  double x_len = 1.0;
  double y_len = 1.0;
  pMeshEnt vertices[(x+1) * (y+1)]; 

  //create the vertices at the assigned positions
  for (unsigned int i = 0; i < (x+1) * (y+1); ++i)
  {
    double points[3];
    points[0] = (i - i/(x+1)*(x+1))*x_len/x;
    points[1] = (i/(x+1))*y_len/y;
    points[2] = 0.0;
    vertices[i] = pumi_mesh_createVtx(mesh, NULL, points);
    double temp[2];
    pumi_node_getCoord(vertices[i], 0, temp);
    std::cout << "vertix " << i << " coordinate is (" << temp[0] << ", " << temp[1] << ")" << std::endl; 
  }

  pMeshEnt vertices_selected[3];
  vertices_selected[0] = vertices[0];
  vertices_selected[1] = vertices[1];
  vertices_selected[2] = vertices[2];
  pMeshEnt e_ = pumi_mesh_createElem(mesh, NULL, PUMI_TRIANGLE, vertices_selected);
  std::cout << "element 0 " << pumi_ment_getID(e_) << " created" << std::endl;

  vertices_selected[0] = vertices[1];
  vertices_selected[1] = vertices[3];
  vertices_selected[2] = vertices[2];
  e_ = pumi_mesh_createElem(mesh, NULL, PUMI_TRIANGLE, vertices_selected);
  std::cout << "element 1 " << pumi_ment_getID(e_) << " created" << std::endl;
  
#if 1

  pMeshEnt vertices4;
  pMeshEnt vertices5;

  double points[3];
  points[0] = 2.0;
  points[1] = 0.0;
  points[2] = 0.0;

  vertices4 = pumi_mesh_createVtx(mesh, NULL, points);

  points[1] = 1.0;
  vertices5 = pumi_mesh_createVtx(mesh, NULL, points);

  pMeshEnt vertices_quad[4];
  vertices_quad[0] = vertices[1];
  vertices_quad[1] = vertices4;
  vertices_quad[2] = vertices5;
  vertices_quad[3] = vertices[3];

  std::cout << "all vertices created" << std::endl;

  e_ = pumi_mesh_createElem(mesh, NULL, PUMI_QUAD, vertices_quad);
  std::cout << "element 2 " << pumi_ment_getID(e_) << " created" << std::endl;

#endif



#if 0
  //create the element (faces) at the assigned positions
  //the order of vertices follow the anti-clockwise order mentioned in the manual
  for (unsigned int i = 0; i < y; ++i)
  { 
    for (unsigned int j = 0; j < x; ++j)
    {
      pMeshEnt vertices_selected[4];
      int k = i * x + j; //index of element
      vertices_selected[0] = vertices[k + i];
      vertices_selected[1] = vertices[k + i + 1];
      vertices_selected[2] = vertices[k + i + 8];
      vertices_selected[3] = vertices[k + i + 7];
      pMeshEnt e = pumi_mesh_createElem(mesh, NULL, PUMI_QUAD, vertices_selected);
      std::cout << "element " << k << " " << pumi_ment_getID(e) << " created" << std::endl;
    }
  }
#endif

  pMeshEnt e;
  pMeshIter it = mesh->begin(2);
  //iterate all the faces and print out the adjacent vertices
  //to verify whether the faces are created correctly
  while ((e = mesh->iterate(it)))
  {
    std::vector<pMeshEnt> v;
    pumi_ment_getAdj(e, 0, v);
    std::cout << "element " << pumi_ment_getID(e) << " has vertices "; 
    for (unsigned int i = 0; i < v.size(); ++i)
      std::cout << pumi_ment_getID(v[i]) << " ";
    std::cout << std::endl;
  }


  std::string ele_type(argv[1]);
  if (ele_type.compare("quadratic") == 0)
      pumi_mesh_setShape(mesh, pumi_shape_getSerendipity());

  Reorder(mesh, ele_type);


  pumi_mesh_freeze(mesh);
  pumi_mesh_verify(mesh);

  char mesh_name[100];
  strcpy(mesh_name, argv[2]);
  strcat(mesh_name, ".smb");
      

  pumi_mesh_write(mesh,mesh_name);
  pumi_mesh_write(mesh,argv[2],"vtk");
  pumi_mesh_delete(mesh);
  pumi_finalize();
  MPI_Finalize();
}

