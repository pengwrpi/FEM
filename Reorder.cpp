#include <list>
#include <algorithm>
#include "FEA.h"

void enqueue_list(std::list<pMeshEnt>& equeue, pMeshEnt& e) {
    //determine if the e already exists
    //if not, push it back
    if (equeue.empty())
        equeue.push_back(e);
    else if (std::find(equeue.begin(), equeue.end(), e) == equeue.end())
        equeue.push_back(e);
}

bool hasNode(pMesh m, pMeshEnt e) {
  return pumi_shape_hasNode(pumi_mesh_getShape(m),pumi_ment_getTopo(e));
}

pMeshEnt Reorder_getstart(pMesh m) {
// find the vertex with least number of edges bounded
  int min_edge; //minimum number of edges that a vertex can have
  int index_v = 0;
  pMeshEnt ent_min_edge; //the entity to be returned
  pMeshIter it = m->begin(0);
  pMeshEnt e;
  while ((e = m->iterate(it))) {
    int num_edge;
    if (index_v == 0) {
      //initialize the min_edge and ent_min_edge
      min_edge = pumi_ment_getNumAdj(e, 1);
      ent_min_edge = e;
    } else {
      num_edge = pumi_ment_getNumAdj(e, 1); 
      if (num_edge < min_edge) {
        min_edge = num_edge;
        ent_min_edge = e;
      }
    }  
    if (min_edge == 2) {
      return ent_min_edge;
    }
    index_v++;
  }
  m->end(it);
  return ent_min_edge;
}  
  

void FEA::Reorder() {

  int labelnode;//label of node
  int labelface = pumi_mesh_getNumEnt(mesh, 2) + 1;//label of faces
  if(ele_type.compare("quadratic") == 0) {
    //if edge midpoint has a node
    labelnode = pumi_mesh_getNumEnt(mesh, 0) + pumi_mesh_getNumEnt(mesh, 1);
  } else {
    labelnode = pumi_mesh_getNumEnt(mesh, 0);
  }

  numbers = pumi_numbering_create(mesh, "my_numbers", NULL); //create an empty numbering for reordering

  pMeshEnt e = Reorder_getstart(mesh); //find the entity with fewest adjacent edges
  std::list<pMeshEnt> equeue; //use STD list to store the entity in queue
  equeue.push_back(e);
  while ( equeue.empty() == false ) {
    pMeshEnt curr_ent = equeue.front();
    equeue.pop_front();//dequeue
    if (!pumi_ment_isNumbered(curr_ent, numbers)) {
      labelnode--;
      pumi_ment_setNumber(curr_ent, numbers, 0, 0, labelnode);
    }
    if (pumi_ment_getDim(curr_ent) == 0) {
      std::vector<pMeshEnt> edges;
      pumi_ment_getAdj(curr_ent, 1, edges);
      std::list<pMeshEnt> To_queue;//initialize the list
      for (unsigned int i = 0; i < edges.size(); ++i) {//loop over the adjacent edges of vertex curr_ent
        pMeshEnt othervertex = pumi_medge_getOtherVtx(edges[i], curr_ent); //oververtex of edge[i]
        if (ele_type.compare("quadratic") == 0) {
          if ( (pumi_ment_isNumbered(othervertex, numbers) || \
                std::find(equeue.begin(), equeue.end(), othervertex) != equeue.end() ) && \
                  !pumi_ment_isNumbered(edges[i], numbers) ) {
            // if the othervertex isn't numbered or in queue and this edges hasn't been numbered
            labelnode--;
            pumi_ment_setNumber(edges[i], numbers, 0, 0, labelnode);
          } else if ( (std::find(equeue.begin(), equeue.end(), othervertex) == equeue.end() || \
                      !pumi_ment_isNumbered(othervertex, numbers) ) &&
                      !pumi_ment_isNumbered(edges[i], numbers) ) {
            // if the othervertex isn't numbered or in queue 
            //equeue.push_back(edges[i]);
            //To_queue.push_back(othervertex);
            enqueue_list(equeue, edges[i]);
            enqueue_list(To_queue, othervertex);
          }
        } else {
          if (!pumi_ment_isNumbered(othervertex, numbers))
            //To_queue.push_back(othervertex);
            enqueue_list(To_queue, othervertex);
        }
      }
        if (!To_queue.empty()) {
          equeue.insert(equeue.end(), To_queue.begin(), To_queue.end());
          To_queue.clear();
        }
    }
  }
  std::cout << "reorder done" << std::endl;
  pumi_mesh_write(mesh, "numbered", "vtk");
}

