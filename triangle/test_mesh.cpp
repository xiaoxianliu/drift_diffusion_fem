#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include "mesh.hpp"
using namespace std;

/* 1. Generate mesh and compute geometric and topological quantities
/* 2. Read mesh(nodes, edges, elements, topologies)
/* 3. Plot mesh 
*/

#define OUTPUT output

int main(){
	MeshData mesh;
	/* Generate mesh by running a script calling "triangle" */
	/*	# q: quality mesh
	/*	# p: read input from a .poly file
	/*	# e: generate edge file
	/*	# A: apply regional attributes
	/*	# a: area constraint
	*/
	string polyfile_name = "example";

	string cmd="/home/xiaoxian/bin/triangle/triangle -qzpeAa0.001 " + polyfile_name + ".poly";
	system(cmd.c_str());

	/* Read in mesh information and compute mesh quantities */
	string node_file_name= polyfile_name + ".1.node";
	string edge_file_name= polyfile_name + ".1.edge";
	string ele_file_name = polyfile_name + ".1.ele";

	ReadNodes(mesh, node_file_name);
	ReadEdges(mesh, edge_file_name);
	ReadElements(mesh, ele_file_name);

	ComputeTopology(mesh);

	ComputeEdgeLengths(mesh);
	ComputeElementAreas(mesh);

	WriteGNUplot(mesh, polyfile_name);


	/****************************** OUTPUT ************************************/
	/* Write output file */
	ofstream output;
	output.open("output.txt");

	OUTPUT << "Nodal information:" << endl;
	for (int i=0; i<mesh.num_nodes; i++)
	{
		OUTPUT << i << "\t" << mesh.nodes[i*2] << "\t" << mesh.nodes[i*2+1] << "\t" << mesh.node_markers[i] << endl;
	}
	OUTPUT << "\n";

	OUTPUT << "Edge information:" << endl;
	for (int i=0; i<mesh.num_edges; i++)
	{
		OUTPUT << i << "\t" << mesh.edges[i*2] << "\t" << mesh.edges[i*2+1] << "\t" << mesh.edge_markers[i] << endl;
	}
	OUTPUT << "\n";

	OUTPUT << "Element information:" << endl;
	for (int i=0; i<mesh.num_elements; i++)
	{
		OUTPUT << i << "\t";
		for (int j=0; j<mesh.num_nodes_per_ele; j++)
			{OUTPUT << mesh.elements[ i*mesh.num_nodes_per_ele + j ] << "\t";	}
		OUTPUT << mesh.element_markers[i] << endl << endl;
	}
	OUTPUT << endl;


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
	OUTPUT << "\n\n";

	/************ Topology ends */

	/*********** Geometric quantities ***********/
	OUTPUT << "Lengths of edges" << endl;
	for (int i=0; i<mesh.num_edges; i++)
		OUTPUT << "Edge "<<i<<":\t" << mesh.edge_lengths[i] << "\n";
	cout << "\n";

	OUTPUT << "Areas of elements" << endl;
	for (int i=0; i<mesh.num_elements; i++)
		OUTPUT << "Element " << i << ":\t" << mesh.ele_areas[i] << "\n";
	cout << "\n\n";



	

	output.close();
	return 0;

}
