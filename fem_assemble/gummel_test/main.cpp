#include <iostream>
#include <vector>
#include <armadillo>
#include "../../triangle/mesh.hpp"
#include "../fem_assemble.hpp"

using namespace std;
using namespace my_mesh;
using namespace linear_fem;

int main(int argc, char* argv[])
{

	// 1. Define mesh
	vector< vector<double> > interface_nodes;
	vector<double> new_node(2);
	new_node[0]=0;	new_node[1]=1;	interface_nodes.push_back(new_node);
	new_node[0]=0;	new_node[1]=0;	interface_nodes.push_back(new_node);

	double max_area=0.001;
	MeshData mesh = generateMesh("gummel_test", interface_nodes, max_area);
	
	return 0;
}
