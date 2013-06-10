#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "../mesh.hpp"


/* 1. Generate mesh and compute geometric and topological quantities. Plot mesh
/* 2. Refine mesh -> new_mesh, and compute its topological and geometric attributes and plot it.
*/

using namespace std;
using namespace MeshNamespace;

int main(){



	/* 1. Form new mesh input */
	vector<vector<double> > interface_nodes;
	vector<int> node_markers;
	vector<vector<int> > segments;
	vector<int> segment_markers;
	vector<vector<double> > regions;
	vector<int> region_markers;

	/* 1.1 Add interface node to the vector of interface nodes */
	vector<double> new_node(2);
	new_node[0] = 0.0;	new_node[1]=1.0;
	interface_nodes.push_back(new_node);
	new_node[0] = 0.0;	new_node[1]=0.0;
	interface_nodes.push_back(new_node);


	newMesh(interface_nodes, node_markers,		// 4 corner vertices are added to "interface_nodes"
		segments, segment_markers,
		regions, region_markers);

	vector<vector<double> > nodes = interface_nodes;

	/* 2. Write new mesh info to a .poly file */
	string filename = "rectangle";


	writePolyfile(	filename,
			nodes, node_markers,			// vertices
			segments, segment_markers,		// segments
			regions, region_markers);		// regional attributes




	/* 3. Generate mesh by running a script calling "triangle" */
	/*	# q: quality mesh
	/*	# p: read input from a .poly file
	/*	# e: generate edge file
	/*	# A: apply regional attributes
	/*	# a: area constraint
	*/
	string cmd = "/home/xiaoxian/bin/triangle/triangle -qzpeAa1 " + filename + ".poly";
	system(cmd.c_str());



	/* 4. Read in mesh information and compute mesh quantities */
	MeshData mesh;

	string node_file_name= filename + ".1.node";
	string edge_file_name= filename + ".1.edge";
	string ele_file_name = filename + ".1.ele";

	ReadNodes(mesh, node_file_name);
	ReadEdges(mesh, edge_file_name);
	ReadElements(mesh, ele_file_name);

	ComputeTopology(mesh);

	ComputeEdgeLengths(mesh);
	ComputeElementAreas(mesh);

	/* 4.1 Get info of interface */
	vector<int> interface_nodes_extracted;
	vector<int> interface_edges_extracted;
	extractInterface(mesh, interface_edges_extracted, interface_nodes_extracted);

	/* 5. Plot mesh */
	gnuplot_mesh(mesh, filename);
	gnuplot_interface(mesh, filename);



	/* 6 Refine mesh */
	/* 6.1 Form marker for refinement*/
	mesh.refinement_markers.resize(mesh.num_elements, 0);
	mesh.refinement_markers[0] = 1;

	/* 6.2 refine mesh */
	refineMesh(mesh, filename, 1);


	return 0;

}
