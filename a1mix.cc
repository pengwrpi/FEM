#include <pumi.h>
#include <apf.h>
#include <PCU.h>
#include <iostream>
#include <Eigen/Dense>
#include "FEA.h"


int adjacency_mesh_search(/* in */ pMesh m,
			  /* in */ int *initial_simplex,
			  /* in */ double *final_position,
                          /* out */ int *final_count,
			  /* out */ int *final_simplex);

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  pGeom geom = pumi_geom_load("model.dmg", "mesh");
  pMesh mesh = pumi_mesh_load(geom, "part.smb", 4);

  pMeshEnt e;
  pMeshIter it = mesh->begin(2);
  while ((e = mesh->iterate(it)))
  {
    int initialid = pumi_ment_getID(e);
    double positions[2];
    positions[0] = 2.0;
    positions[1] = 0.0;
    int finalid;
    int count_;
    adjacency_mesh_search(mesh, &initialid, positions, &count_, &finalid);

    std::cout << "idout = " << finalid << "for element " << initialid << std::endl;
    break;
  }
  mesh->end(it);

  pumi_mesh_write(mesh,"outMixed","vtk");
  pumi_mesh_delete(mesh);
  pumi_finalize();
  MPI_Finalize();
}


int adjacency_mesh_search(/* in */ pMesh m,
			  /* in */ int *initial_simplex,
			  /* in */ double *final_position,
			  /* out */ int *final_count,
			  /* out */ int *final_simplex)
//*******************************************************
{
  bool located = false;
  pMeshEnt e = NULL;
  pMeshEnt simplex = NULL;
  apf::Adjacent adjacent;
  int edge_curr_index;  
  int edge_prev_index = -1;
  int simplex_dim = pumi_mesh_getDim(m); 
  int vertex_dim = 0, edge_dim = 1;
  int bridge_dim = simplex_dim - 1;
  int count = 0, max_count = pumi_mesh_getNumEnt(m, simplex_dim); 
  double tol = 1e-15;
  
  simplex = pumi_mesh_findEnt(m, simplex_dim, *initial_simplex);
  while (not located) {
    int simplex_index = pumi_ment_getID(simplex);//doubt
    std::vector<pMeshEnt> vertices;
    Vector3 v_coord[3];
    apf::Matrix<3, 3> A;
    pumi_ment_getAdj(simplex, vertex_dim, vertices);
    int nv = vertices.size();
    for (int j = 0; j < nv; ++j) {
      pumi_node_getCoordVector(vertices[j], 0, v_coord[j]);
      // Note: second argument is 0 for linear meshes
    }
    // Compute (linear) barycentric coordinates of final position
    // with respect to current simplex
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 3; ++k)
	A[j][k] = v_coord[k][j];
    for (int j = 0; j < 3; ++j)
      A[2][j] = 1;
    apf::Matrix<3, 3> Ainv = apf::invert(A);
    apf::Vector3 b_coords; // = Ainv * final_position
    for (int j = 0; j < 3; ++j) {
      b_coords[j] = (Ainv[j][0] * final_position[0] +
		     Ainv[j][1] * final_position[1] +
		     Ainv[j][2]);
    }
    // If all positive for current simplex, exit.
    if (((b_coords[0] >= -tol) && (b_coords[0] <= (1 + tol))) &&
	((b_coords[1] >= -tol) && (b_coords[1] <= (1 + tol))) &&
	((b_coords[2] >= -tol) && (b_coords[2] <= (1 + tol)))) {
      located = true;
      *final_simplex = simplex_index;
      break;
    }
    // Otherwise, check which coordinates are negative
    //bool b_negative[3] = {false, false, false};
    int bneg_index[2] = {0, 0};
    int bneg_count = 0;
    for (int j = 0; j < 3; ++j)
      if (b_coords[j] < 0) {
	//b_negative[j] = true;
	bneg_index[bneg_count] = j;
	++bneg_count;
      }
    // Obtain the index of most negative coordinate
    //assert(bneg_count > 0 && bneg_count < 3);

    // Ensure bneg_index[0] is the index of vertex whose corresponding
    // barycentric coordinate is most negative; ties automatically resolved
    // as a result.
    if (bneg_count == 2)
      if (fabs(b_coords[bneg_index[0]]) <
	  fabs(b_coords[bneg_index[1]])) {
	int tmp_index = bneg_index[0];
	bneg_index[0] = bneg_index[1];
	bneg_index[1] = tmp_index;
      }

    // Determine edge opposite to this vertex
    pMeshEnt edge_vertices[2];
    int edge_count = 0;
    int edge_type = 1;    // Mesh entity type; see APF documentation.
    for (int j = 0; j < 3; ++j)
      if (j != bneg_index[0]) {
	edge_vertices[edge_count] = vertices[j];
	++edge_count;
      }
    
    // Find neighboring simplex sharing this edge
    e = apf::findElement(m, edge_type, edge_vertices);
    edge_curr_index = pumi_ment_getID(e);

    // If current edge choice is same as previous edge, pick edge
    // opposite to second least (actual, not absolute, valued) barycentric
    // coordinate.
    if (edge_curr_index == edge_prev_index) {
      edge_count = 0;
      for (int j = 0; j < 3; ++j)
	if (j != bneg_index[1]) {
	  edge_vertices[edge_count] = vertices[j];
	  ++edge_count;
	}
      e = apf::findElement(m, edge_type, edge_vertices);
      edge_curr_index = pumi_ment_getID(e);
    }     

    apf::getBridgeAdjacent(m, e,
			   bridge_dim, simplex_dim,
			   adjacent);
    if (adjacent.getSize() == 2) {
      for (size_t j = 0; j < adjacent.getSize(); ++j) {
	int new_simplex_index = pumi_ment_getID(adjacent[j]);
	if (new_simplex_index != simplex_index)
	  simplex = adjacent[j];
      }
    }
    else {
      apf::Downward edges;
      int ne = m->getDownward(simplex, edge_dim, edges);
      for (int j = 0; j < ne; ++j) {
	int edge_tmp_index = pumi_ment_getID(edges[j]);
	if ((edge_tmp_index != edge_curr_index) &&
	    (edge_tmp_index != edge_prev_index)) {
	  e = edges[j];
	  break;
	}
      }
      edge_curr_index = pumi_ment_getID(e);    
      apf::getBridgeAdjacent(m, e,
			     bridge_dim, simplex_dim,
			     adjacent);
      for (size_t j = 0; j < adjacent.getSize(); ++j) {
	int new_simplex_index = pumi_ment_getID(adjacent[j]);
	if (new_simplex_index != simplex_index)
	  simplex = adjacent[j];
      }
    }
    
    // Keep track of edge via which we entered the current simplex
    edge_prev_index = edge_curr_index;
    ++count;

    if (count == max_count){
      // std::cout << "\nError: hit maximum number of simplices to look for.";
      *final_simplex = -2;
      *final_count = count;
      return 1;

    }
  }
  *final_count = count;
  return 0;
}
