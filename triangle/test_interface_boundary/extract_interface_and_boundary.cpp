#include <iostream>
#include <string>
#include <vector>
#include "../mesh.hpp"

int main()
{
using namespace std;
using namespace my_mesh;

	// 1. Generate mesh
	string filename = "test_extract_boundary_and_interface";

	vector< vector<double> > interface_nodes;
	vector<double> new_node(2);
	new_node[0] = 0.0;	new_node[1]=1.0;
	interface_nodes.push_back(new_node);
	new_node[0] = 0.0;	new_node[1]=0.0;
	interface_nodes.push_back(new_node);

	double max_area = 0.1;
	MeshData mesh;
	mesh = generateMesh(filename, interface_nodes, max_area );

	// 2. Output boundary and interface information
	// 2.1 boundary "1"
	for (int i=0; i<mesh.boundary1_nodes.size(); i++)
	{	int node_index = mesh.boundary1_nodes[i];
		cout << i << ". Node " << node_index << "\t" 
			<< mesh.nodes[node_index][0] << "\t" << mesh.nodes[node_index][1] << "\n";
	}
	for (int i=0; i<mesh.boundary1_edges.size(); i++)
	{	int edge_index = mesh.boundary1_edges[i];
		cout << i << ". Edge " << edge_index << "\t"
			<< mesh.edges[edge_index][0] << "\t" << mesh.edges[edge_index][1] << "\n";
	}
	cout << "\n\n";

	// 2.2 boundary "3"
	for (int i=0; i<mesh.boundary3_nodes.size(); i++)
	{	int node_index = mesh.boundary3_nodes[i];
		cout << i << ". Node " << node_index << "\t" 
			<< mesh.nodes[node_index][0] << "\t" << mesh.nodes[node_index][1] << "\n";
	}
	for (int i=0; i<mesh.boundary3_edges.size(); i++)
	{	int edge_index = mesh.boundary3_edges[i];
		cout << i << ". Edge " << edge_index << "\t"
			<< mesh.edges[edge_index][0] << "\t" << mesh.edges[edge_index][1] << "\n";
	}
	cout << "\n\n";

	// 2.3 interface "5"
	// 2.1 boundary "1"
	for (int i=0; i<mesh.interface_nodes.size(); i++)
	{	int node_index = mesh.interface_nodes[i];
		cout << i << ". Node " << node_index << "\t" 
			<< mesh.nodes[node_index][0] << "\t" << mesh.nodes[node_index][1] << "\n";
	}
	for (int i=0; i<mesh.interface_edges.size(); i++)
	{	int edge_index = mesh.interface_edges[i];
		cout << i << ". Edge " << edge_index << "\t"
			<< mesh.edges[edge_index][0] << "\t" << mesh.edges[edge_index][1] << "\n";
	}


return 0;
}
