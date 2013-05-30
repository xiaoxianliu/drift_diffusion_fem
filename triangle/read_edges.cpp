#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "mesh.hpp"

using namespace std;

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
	int num_edges, num_marker_per_edge;
	linestream >> num_edges >> num_marker_per_edge;
	if (num_marker_per_edge!=0 and num_marker_per_edge!=1)
		{cout << "Number of edge markers is " << num_marker_per_edge << ".\n";
		 cout << "It has to be either 0 or 1. \n";
		 exit(1);}

	mesh.num_edges = num_edges;
	mesh.num_marker_per_edge = num_marker_per_edge;


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
		if (num_marker_per_edge==1)
			linestream >> marker;

		/* Finally, set cooresponding values of "mesh" to the temporary variables above*/
		mesh.edges.push_back(edge);
		mesh.edge_markers.push_back(marker);

	}
	input_file.close();

	/* Final check on the total number of edges */
	if (mesh.num_edges != mesh.edges.size())
		{ cout << "\"mesh.num_edges\" is " << mesh.num_edges << "\n";
		  cout << "\"mesh.edges\" is a cpp vector of size " << mesh.edges.size() << "\n";
		  cout << "They have to be equal! Exit... \n";
		  exit(1);
		}

	return 0;
}
