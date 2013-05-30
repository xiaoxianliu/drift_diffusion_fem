#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>

#include "mesh.hpp"

/* Read in information of edges */
/* Input: (1) MeshData &mesh;
	  (2) string edge_file_name.
*/

int ReadEdges(MeshData &mesh, string edge_file_name){
	ifstream input_file;
	input_file.open(edge_file_name.c_str());
	if (!input_file.is_open())
		{cout<< "Faile to open edge file " << edge_file_name << endl;	exit(1);}

	/* Read in: # of edges; # of edge marker (0 or 1 by default). */
	string line;
	getline(input_file, line);
	stringstream linestream(line);
	int num_edges, num_markers;
	linestream >> num_edges >> num_markers;
	if (num_markers!=0 and num_markers!=1)
		{cout << "Number of edge markers is " << num_markers << ".\n";
		 cout << "It has to be either 0 or 1. \n";
		 exit(1);}


	/* allocate memory for relavent element of "mesh" */
	mesh.num_edges = num_edges;


	/* Read in edge information: nodes belong to each edge; edge marker */
	while (input_file.good())
	{
		getline(input_file, line);

		/* Skip unimportant lines */
		if (line.length()==0)
			{continue;	}				// blank lines
		if (line[0]=='#')
			{continue;	}				// commentary lines

		/* Read in edge information to temporary variables */
		stringstream linestream(line);
		int index;
		vector<int> edge(2);
		int marker=0;

		linestream >> index >> edge[0] >> edge[1];
		if (num_markers==1)
			linestream >> marker;

		/* Finally, set cooresponding values of "mesh" to the temporary variables above*/
		mesh.edges.push_back(edge);
		mesh.edge_markers.push_back(marker);

	}
	input_file.close();
	return 0;
}
