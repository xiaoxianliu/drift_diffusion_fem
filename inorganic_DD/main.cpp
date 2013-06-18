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
	int initial_num_interface_nodes = 401;
	std::vector< std::vector<double> > interface_nodes;
	std::vector<double> new_node(2);
	for (int i=0; i<initial_num_interface_nodes; i++)
	{	new_node[0] = 0;
		new_node[1] = 1 - static_cast<double>(i)/(initial_num_interface_nodes - 1);
		interface_nodes.push_back(new_node);
	}

	my_mesh::MeshData mesh;	
	double max_area = 0.003;
	mesh = my_mesh::generateMesh(filename, interface_nodes, max_area);

	/*** 2. Gummel's map *********/
	arma::vec psi, n, p;
	gummelIteration(mesh, psi, n, p);

	plot_ArmaVec(mesh, psi, "psi");
	plot_ArmaVec(mesh, n, "n");
	plot_ArmaVec(mesh, p, "p");

	

	return 0;
}
