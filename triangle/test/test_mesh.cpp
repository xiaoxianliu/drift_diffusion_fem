#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "../mesh.hpp"
using namespace std;

/* 1. Generate mesh and compute geometric and topological quantities
/* 2. Read mesh(nodes, edges, elements, topologies)
/* 3. Plot mesh 
*/

#define OUTPUT output

int main()
{
using namespace MeshNamespace;

	/* Form new interface */
	vector<vector<double> > interface_nodes;

	vector<double> new_node(2);
	new_node[0] = 0.0;	new_node[1]=1.0;
	interface_nodes.push_back(new_node);
	new_node[0] = 0.0;	new_node[1]=0.0;
	interface_nodes.push_back(new_node);

	/* Generate new mesh */
	MeshData mesh;
	double max_area = 0.001;
	bool is_to_plot=false;
	mesh = generateMesh("test", interface_nodes, max_area);



	/****************************** OUTPUT ************************************/
	/* Write output file */
	ofstream output;
	output.open("test_output.txt");

	OUTPUT << "Nodal information:" << endl;
	for (int i=0; i<mesh.num_nodes; i++)
	{
		OUTPUT << i << "\t" << mesh.nodes[i][0] << "\t" << mesh.nodes[i][1] << "\t" << mesh.node_markers[i] << endl;
	}
	OUTPUT << "\n";

	OUTPUT << "Edge information:" << endl;
	for (int i=0; i<mesh.num_edges; i++)
	{
		OUTPUT << i << "\t" << mesh.edges[i][0] << "\t" << mesh.edges[i][1] << "\t" << mesh.edge_markers[i] << endl;
	}
	OUTPUT << "\n";

	OUTPUT << "Element information:" << endl;
	for (int i=0; i<mesh.num_elements; i++)
	{
		OUTPUT << i << "\t";
		for (int j=0; j<mesh.num_nodes_per_ele; j++)
			{OUTPUT << mesh.elements[i][j] << "\t";	}
		OUTPUT << mesh.element_markers[i] << endl << endl;
	}
	OUTPUT << "\n";


	/************* Topology Start */
	/* 0-to-(0,1,2) topology */
	OUTPUT << "Topology 0 to 0" << endl;
	for (int i=0; i<mesh.num_nodes; i++)
	{	OUTPUT << "Node "<<i<<":\t";
		for (int j=0; j<mesh.topology0to0[i].size(); j++)
			OUTPUT << mesh.topology0to0[i][j]<<"\t";
		OUTPUT << endl;
	}

	OUTPUT << "Topology 0 to 1" << endl;
	for (int i=0; i<mesh.num_nodes; i++)
	{	OUTPUT << "Node "<<i<<":\t";
		for (int j=0; j<mesh.topology0to1[i].size(); j++)
			OUTPUT << mesh.topology0to1[i][j]<<"\t";
		OUTPUT << endl;
	}

	OUTPUT << "Topology 0 to 2" << endl;
	for (int i=0; i<mesh.num_nodes; i++)
	{	OUTPUT << "Node "<<i<<":\t";
		for (int j=0; j<mesh.topology0to2[i].size(); j++)
			OUTPUT << mesh.topology0to2[i][j]<<"\t";
		OUTPUT << endl;
	}
	OUTPUT << "\n\n";

	/* 1-to-(0,1,2) topology */
	OUTPUT << "Topology 1 to 0" << endl;
	for (int i=0; i<mesh.num_edges; i++)
	{	OUTPUT << "Edge " << i << ":\t";
		for (int j=0; j<mesh.topology1to0[i].size(); j++)
		{	OUTPUT << mesh.topology1to0[i][j]<<"\t";
		}
		OUTPUT << endl;
	}

	OUTPUT << "Topology 1 to 1" << endl;
	for (int i=0; i<mesh.num_edges; i++)
	{	OUTPUT << "Edge " << i << ":\t";
		for (int j=0; j<mesh.topology1to1[i].size(); j++)
		{	OUTPUT << mesh.topology1to1[i][j] <<"\t";
		}
		OUTPUT << endl;
	}

	OUTPUT << "Topology 1 to 2" << endl;
	for (int i=0; i<mesh.num_edges; i++)
	{	OUTPUT << "Edge " << i << ":\t";
		for (int j=0; j<mesh.topology1to2[i].size(); j++)
		{	OUTPUT << mesh.topology1to2[i][j] << "\t";
		}
		OUTPUT << endl;
	}
	OUTPUT << "\n\n";


	/* begin 2-to-(0,1,2) topology */
	OUTPUT << "Topology 2 to 0" << endl;
	for (int i=0; i<mesh.num_elements; i++)
	{	OUTPUT << "Element " << i << ":\t";
		for (int j=0; j<mesh.topology2to0[i].size(); j++)
			OUTPUT << mesh.topology2to0[i][j] << "\t";
		OUTPUT << endl;
	}
	
	OUTPUT << "Topology 2 to 1" << endl;
	for (int i=0; i<mesh.num_elements; i++)
	{	OUTPUT << "Element " << i << ":\t";
		for (int j=0; j<mesh.topology2to1[i].size(); j++)
			OUTPUT << mesh.topology2to1[i][j] << "\t";
		OUTPUT << endl;
	}

	OUTPUT << "Topology 2 to 2" << endl;
	for (int i=0; i<mesh.num_elements; i++)
	{	OUTPUT << "Element " << i << ":\t";
		for (int j=0; j<mesh.topology2to2[i].size(); j++)
			OUTPUT << mesh.topology2to2[i][j] << "\t";
		OUTPUT << endl;
	}
	OUTPUT << "\n\n\n";

	/************ Topology ends */

	/*********** Geometric quantities ***********/
	OUTPUT << "Lengths of edges" << endl;
	for (int i=0; i<mesh.num_edges; i++)
		OUTPUT << "Edge "<<i<<":\t" << mesh.edge_lengths[i] << "\n";
	OUTPUT << "\n\n";

	OUTPUT << "Areas of elements" << "\n";
	for (int i=0; i<mesh.num_elements; i++)
		OUTPUT << "Element " << i << ":\t" << mesh.ele_areas[i] << "\n";
	OUTPUT << "\n\n";


	/*********** Submesh begins *******************/
	/* Interface */
	OUTPUT << "Interface nodes are:\n";
	for (int i=0; i<mesh.interface_nodes.size(); i++)
	{	int node_index = mesh.interface_nodes[i];
		OUTPUT << node_index << "\t" << mesh.nodes[node_index][0] << "\t" << mesh.nodes[node_index][1] << "\n";
	}
	OUTPUT << "\n\n";

	OUTPUT << "Interface edges are:\n";
	for (int i=0; i<mesh.interface_edges.size(); i++)
	{	int edge_index = mesh.interface_edges[i];
		OUTPUT << edge_index << "\t" << mesh.edges[edge_index][0] << "\t" << mesh.edges[edge_index][1] << "\n";
	}

	/************ Submesh ends *******************/

	output.close();


	return 0;

}
