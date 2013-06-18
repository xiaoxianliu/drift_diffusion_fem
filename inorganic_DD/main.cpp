#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include "../triangle/mesh.hpp"
#include "inorganicDD.hpp"


int main (int argc, char* argv[])
{
	std::string filename = "dd_model";

	/*** 1. Mesh ******************/
	std::vector< std::vector<double> > interface_nodes;
	std::vector<double> new_node(2);
	new_node[0]=0;	new_node[1]=1;	interface_nodes.push_back(new_node);
	new_node[0]=0;	new_node[1]=0;	interface_nodes.push_back(new_node);

	my_mesh::MeshData mesh;	
	double max_area = 1;
	mesh = my_mesh::generateMesh(filename, interface_nodes, max_area);

	/*** 2. Gummel's map *********/
	arma::vec psi, n, p;
	gummelIteration(mesh, psi, n, p);

	

	return 0;
}
