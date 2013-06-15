#include <iostream>
#include <cstdlib>
#include <vector>

#include "mesh.hpp"

int MeshNamespace::computeBarryPoints( MeshNamespace::MeshData &mesh)
{
	// Safety checking
	if (mesh.num_elements<=0 || mesh.num_nodes <=0 || mesh.topology2to0.size() <=0)
	{	std::cout << "Number of elements is " << mesh.num_elements << "\n";
		std::cout << "Number of nodes is " << mesh.num_nodes << "\n";
		std::cout << "Length of \"topology2to0\" is " << mesh.topology2to0.size() << "\n";
		std::cout << "Have to construct mesh first! \n";
		std::exit(1);
	}

	if (mesh.topology2to0.size() != mesh.num_elements )
	{	std::cout << "\"Topology2to0\" and number of elements should be equal!\n";
		std::exit(1);
	}

	// Reset STL vector
	mesh.barry_points.clear();

	// Compute coordinates for each barry point
	for (int t=0; t< mesh.topology2to0.size(); t++)
	{	std::vector<double> barry_pt;

		int v0 = mesh.topology2to0[t][0], v1 = mesh.topology2to0[t][1], v2=mesh.topology2to0[t][2];

		double barry_x, barry_y;
		barry_x = ( mesh.nodes[v0][0] + mesh.nodes[v1][0] + mesh.nodes[v2][0] )/3.0;
		barry_y = ( mesh.nodes[v0][1] + mesh.nodes[v1][1] + mesh.nodes[v2][1] )/3.0;

		barry_pt.push_back(barry_x);
		barry_pt.push_back(barry_y);

		mesh.barry_points.push_back(barry_pt);
	}

	return 0;
}
