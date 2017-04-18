#include <PCU.h>
#include <pumi.h>
#include <iostream>
#include "createmesh.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  pGeom geom = pumi_geom_load(NULL, "null");
  pMesh mesh = pumi_mesh_create(geom, 2);

  //for this specific problem, we know there are 5*7 = 35 vertices
  int x = 100;
  int y = 2;
  pMeshEnt vertices[x * y]; 

  //create the vertices at the assigned positions
  for (unsigned int i = 0; i < x * y; ++i)
  {
    double points[3];
    points[0] = (i - i/x*x)*1.0;
    points[1] = (i/x)*1.0;
    points[2] = 0.0;
    vertices[i] = pumi_mesh_createVtx(mesh, NULL, points);
    double temp[3];
    pumi_node_getCoord(vertices[i], 0, temp);
    std::cout << "vertix " << i << " coordinate is (" << temp[0] << ", " << temp[1] << ")" << std::endl; 
  }

  //create the element (faces) at the assigned positions
  //the order of vertices follow the anti-clockwise order mentioned in the manual
  for (unsigned int i = 0; i < y - 1; ++i)
  { 
    for (unsigned int j = 0; j < x - 1; ++j)
    {
      pMeshEnt vertices_selected[4];
      int k = i * (x - 1) + j; //index of element
      vertices_selected[0] = vertices[k + i];
      vertices_selected[1] = vertices[k + i + 1];
      vertices_selected[2] = vertices[k + i + x + 1];
      vertices_selected[3] = vertices[k + i + x];
      pMeshEnt e = pumi_mesh_createElem(mesh, NULL, PUMI_QUAD, vertices_selected);
      std::cout << "element " << k << " " << pumi_ment_getID(e) << " created" << std::endl;
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

  Reorder(mesh, "linear");

      

  pumi_mesh_write(mesh,"outQuad.smb");
  pumi_mesh_write(mesh,"outQuad","vtk");
  pumi_mesh_delete(mesh);
  pumi_finalize();
  MPI_Finalize();
}

