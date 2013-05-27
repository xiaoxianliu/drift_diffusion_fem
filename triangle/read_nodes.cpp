#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>

#include "mesh.hpp"

using namespace std;

/* Function that modify the node information of given MeshData object "mesh" */
/* Input: (1) MeshData &mesh; (2) full name of ".node" file */
int ReadNodes(MeshData &mesh, string node_file_name){
	ifstream input_file;
	input_file.open(node_file_name.c_str());
	if (!input_file.is_open())
		{cout << "Faile to open file: " << node_file_name << endl;	exit(1);}

	/* Read first line of ".node" file */
	string line;
	getline(input_file, line);
	stringstream linestream(line);
	/* Read in: # of nodes; dimension; # of attributes; # of markers */
	int dimension, num_attributes, num_markers;
	linestream >> mesh.num_nodes >> dimension >> num_attributes >> num_markers;
	if (dimension!=2)
		{cout << "Dimension of mesh has to be 2" << endl; exit(1);}

	/* Read in the nodal information of mesh: node index, node coordinates, node attributes (if any), node marker (if any) */
	size_t index;
	mesh.nodes = new double[mesh.num_nodes*2];
	mesh.node_markers = new int[mesh.num_nodes];

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
		double x, y;
		double attributes[num_attributes];
		int marker=0;
		linestream >> index >> x >> y;				// read index and coordinates of current node
		for (int i=0; i<num_attributes; i++)			// skip all the attributes, if any
			{linestream >> attributes[i];	}
		if (num_markers==1)					// read marker of node
			{linestream >> marker;}

		/* Finally, set the corresponding items in "mesh" to the temporary variables above*/
		mesh.nodes[index*2] = x;	mesh.nodes[index*2 + 1] = y;
		mesh.node_markers[index] = marker;
	}

	input_file.close();	
	return 0;
}
