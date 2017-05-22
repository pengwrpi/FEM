#include <PCU.h>
#include <pumi.h>
#include <iostream>
#include <cstdlib>
#include "createmesh.h"

int main(int argc, char** argv)
{
  if (argc != 7)
  {
      perror("wrong number of argv");
      exit(1);
  }
  MPI_Init(&argc,&argv);
  pumi_start();
  pGeom geom = pumi_geom_load(NULL, "null");
  pMesh mesh = pumi_mesh_create(geom, 2);

  //for this specific problem, we know there are 5*7 = 35 vertices
  int x = atoi(argv[1]);
  int y = atoi(argv[2]);
  double len_x = atof(argv[3]);
  double len_y = atof(argv[4]);
  pMeshEnt vertices[x * y]; 
  double dx = len_x/(x-1);
  double dy = len_y/(y-1);

  //create the vertices at the assigned positions
  for (unsigned int i = 0; i < x * y; ++i)
  {
    double points[3];
    points[0] = (i - i/x*x)*dx;
    points[1] = (i/x)*dy;
    points[2] = 0.0;
    vertices[i] = pumi_mesh_createVtx(mesh, NULL, points);
    double temp[3];
    pumi_node_getCoord(vertices[i], 0, temp);
    std::cout << "vertix " << i << " coordinate is (" << temp[0] << ", " << temp[1] << ")" << std::endl; 
  }

  //create the element (faces) at the assigned positions
  //the order of vertices follow the anti-clockwise order mentioned in the manual
  int k = 0;
  for (unsigned int i = 0; i < y - 1; ++i)
  { 
    for (unsigned int j = 0; j < x - 1; ++j)
    {
      pMeshEnt vertices_selected[3];
      int start = j + i * x; //index of first element
      vertices_selected[0] = vertices[start];  
      vertices_selected[1] = vertices[start + 1];
      vertices_selected[2] = vertices[start + x];
      pMeshEnt e = pumi_mesh_createElem(mesh, NULL, PUMI_TRIANGLE, vertices_selected);
      std::cout << "element " << k << " " << pumi_ment_getID(e) << " created" << std::endl;
      k++;

      start++; //index of second element
      vertices_selected[0] = vertices[start]; 
      vertices_selected[1] = vertices[start + x];
      vertices_selected[2] = vertices[start + x - 1];
      e = pumi_mesh_createElem(mesh, NULL, PUMI_TRIANGLE, vertices_selected);
      std::cout << "element " << k << " " << pumi_ment_getID(e) << " created" << std::endl;
      k++;
    }
  }

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


  pumi_mesh_freeze(mesh);
  pumi_mesh_verify(mesh);

  if (strcmp(argv[5], "quadratic") == 0)
      pumi_mesh_setShape(mesh, pumi_shape_getSerendipity());

  Reorder(mesh, argv[5]);

      
  char mesh_name[100];
  strcpy(mesh_name, argv[6]);
  strcat(mesh_name, ".smb");

  pumi_mesh_write(mesh,mesh_name);
  pumi_mesh_write(mesh,argv[6],"vtk");
  pumi_mesh_delete(mesh);
  pumi_finalize();
  MPI_Finalize();
}

