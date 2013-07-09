#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "mesh.hpp"

#define PI 3.14159265359

namespace MeshNamespace{

/********************************************************************************************************************************/
// Cosine interface
MeshNamespace::MeshData generateMesh_cosine_interface (	std::string filename,
							double x_offset,
							double x_amplitude,
							int num_periods,
							int num_interface_nodes,
							double max_area){
	// 1. Mesh file name
	std::string mesh_filename = filename + "_cos_interface";

	// 2. Interface nodes
	std::vector< std::vector<double> > interface_nodes;
	std::vector<double> node(2);		// coordinate information of each added node

	// 2.1 first end node on interface
	node[1] = 0.0;	node[0] = x_offset;
	interface_nodes.push_back(node);
	// 2.2 intermediate nodes on interface
	for (int i=1; i<num_interface_nodes-1; i++)
	{	double x,y;
		y = (1.0/(num_interface_nodes-1)) * i; 
		x = x_offset + x_amplitude * (1.0 - std::cos(2*num_periods*PI*y)); 

		node[0] = x;
		node[1] = y;
		interface_nodes.push_back(node);
	}
	// 2.3 Last end node on interface
	node[0] = x_offset;	node[1] = 1.0;
	interface_nodes.push_back(node);

	// 3. Finally, generate mesh
	MeshNamespace::MeshData mesh = MeshNamespace::generateMesh(mesh_filename, interface_nodes, max_area, true);

return mesh;
}



/********************************************************************************************************************************/
// "Great wall"/"smart" interface
MeshNamespace::MeshData generateMesh_great_wall (	std::string filename,
							double x_offset,
							double x_amplitude,
							int num_bumps,
							double max_area){

	// 1. Mesh file name
	std::string mesh_filename = filename + "_great_wall";

	// 2. Interface nodes
	int num_interface_nodes = 2 + num_bumps*4;
	std::vector< std::vector<double> > interface_nodes;
	std::vector<double> node(2);		// coordinate information of each added node

	// 2.1 first end node on interface
	node[1] = 0.0;	node[0] = x_offset;
	interface_nodes.push_back(node);
	// 2.2 intermediate nodes on interface
	for (int i=0; i<num_bumps; i++)
	{
		node[0] = x_offset;	node[1] = (1.0/(2.0*num_bumps + 1.0)) * (2*i+1);
		interface_nodes.push_back(node);
		node[0] = x_offset + x_amplitude;
		interface_nodes.push_back(node);

		node[0] = x_offset + x_amplitude;	node[1] = (1.0/(2.0*num_bumps + 1.0)) * (2*i+2);
		interface_nodes.push_back(node);
		node[0] = x_offset;
		interface_nodes.push_back(node);
	}
	// 2.3 Last end node on interface
	node[0] = x_offset;	node[1] = 1.0;
	interface_nodes.push_back(node);

	// 3. Finally, generate mesh
	MeshNamespace::MeshData mesh = MeshNamespace::generateMesh(mesh_filename, interface_nodes, max_area, true);

return mesh;
}

}
