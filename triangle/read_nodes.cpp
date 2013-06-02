#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "mesh.hpp"

using namespace std;
using namespace MeshNamespace;

/* Function that modify the node information of given MeshData object "mesh" */
/* Input: (1) MeshData &mesh; (2) full name of ".node" file */
int MeshNamespace::ReadNodes(MeshData &mesh, string node_file_name){
	ifstream input_file;
	input_file.open(node_file_name.c_str());
	if (!input_file.is_open())
		{cout << "Fail to open file: " << node_file_name << endl;	exit(1);}

	/* 1. Read first line of ".node" file *: general information of all nodes */
	string line;
	getline(input_file, line);
	stringstream linestream(line);
	/* Read in: # of nodes; dimension; # of attributes; # of markers */
	int num_nodes, dimension, num_attributes_per_node, num_marker_per_node;
	linestream >> num_nodes >> dimension >> num_attributes_per_node >> num_marker_per_node;
	if (dimension!=2)
		{cout << "Dimension of mesh has to be 2" << endl; exit(1);}
	mesh.num_nodes = num_nodes;
	mesh.num_attributes_per_node = num_attributes_per_node;
	mesh.num_marker_per_node = num_marker_per_node;



	/* 2. Read in the nodal information of mesh: node index, node coordinates, node attributes (if any), node marker (if any) */
	while (input_file.good())
	{
		/* Skip unimportant lines */
		getline(input_file, line);
		if (line.length()==0)
			{continue;}					//blank line, skip
		if (line[0]=='#')
			{continue;}					//commentary line, skip

		/* Read in nodal information to temporary variables */
		stringstream linestream(line);

		size_t index;
		vector<double> node(2);
		vector<double> attributes(num_attributes_per_node);
		int marker=0;

		linestream >> index >> node[0] >> node[1];		// read index and coordinates of current node
		for (int i=0; i<num_attributes_per_node; i++)			// skip all the attributes, if any
			{linestream >> attributes[i];	}
		if (num_marker_per_node==1)					// read marker of node
			{linestream >> marker;}

		/* Finally, set the corresponding items in "mesh" to the temporary variables above*/
		mesh.nodes.push_back(node);
		mesh.node_attributes.push_back(attributes);
		mesh.node_markers.push_back(marker);
	}

	input_file.close();	

	/* Final check on number of nodes */
	if (mesh.nodes.size() != mesh.num_nodes)
		{ cout << "\"mesh.num_nodes\" is " << mesh.num_nodes << "\n";
		  cout << "size of cpp vector \"mesh.nodes\" is " << mesh.nodes.size() << "\n";
		  cout << "Not equal! Exit...\n";
		  exit(1);
		}

	return 0;
}
