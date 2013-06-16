#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "../mesh.hpp"


/* 1. Generate mesh and compute geometric and topological quantities. Plot mesh
/* 2. Refine mesh -> new_mesh, and compute its topological and geometric attributes and plot it.
*/

using namespace std;
using namespace MeshNamespace;

int main(){

	string filename = "test_refine";

	/* 1. Form new mesh input */
	vector<vector<double> > interface_nodes;
	vector<double> new_node(2);
	new_node[0] = 0.0;	new_node[1]=1.0;
	interface_nodes.push_back(new_node);
	new_node[0] = 0.0;	new_node[1]=0.0;
	interface_nodes.push_back(new_node);

	MeshData mesh = generateMesh(filename, interface_nodes);



	/* 6 Refine mesh */
	/* 6.1 Form marker for refinement*/
	mesh.refinement_markers.resize(mesh.num_elements, 0);
	mesh.refinement_markers[0] = 1;

	/* 6.2 refine mesh */
	refineMesh(mesh, filename, 1);


	return 0;

}
