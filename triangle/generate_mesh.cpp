#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "mesh.hpp"


namespace MeshNamespace
{

MeshData generateMesh(	std::string filename, 
			const std::vector< std::vector<double> > &interface_nodes, 
			double max_area,
			bool is_to_plot)
{
using namespace MeshNamespace;

	MeshData mesh;

	// 1. Initialize inputs for mesh
	std::vector< std::vector<double> > nodes = interface_nodes;	// copy "interface_nodes" to "nodes"
	std::vector<int> node_markers;
	std::vector< std::vector<int> > segments;
	std::vector<int> segment_markers;
	std::vector< std::vector<double> > regions;
	std::vector<int> region_markers;
	initializeMesh (nodes, node_markers, segments, segment_markers, regions, region_markers);
									// Initial boundary nodes are added to "nodes" after this line

	// 2. Generate mesh by "triangle"
	writePolyfile(	filename,
			nodes, node_markers,			// vertices
			segments, segment_markers,		// segments
			regions, region_markers);		// regional attributes


	/* 3. Generate mesh by calling "triangle" */
	/*	# q: quality mesh
	/*	# p: read input from a .poly file
	/*	# e: generate edge file
	/*	# A: apply regional attributes
	/*	# a: area constraint
	*/

	std::stringstream cmd_stream;
	cmd_stream << "/home/xiaoxian/bin/triangle/triangle"
			<< " -qzpeAa" << max_area 
			<< " " + filename + ".poly";
	std::string cmd = cmd_stream.str();
	system(cmd.c_str());



	/* 4. Read in mesh information and compute mesh quantities */
	std::string node_file_name= filename + ".1.node";
	std::string edge_file_name= filename + ".1.edge";
	std::string ele_file_name = filename + ".1.ele";

	ReadNodes(mesh, node_file_name);
	ReadEdges(mesh, edge_file_name);
	ReadElements(mesh, ele_file_name);

	ComputeTopology(mesh);					// mesh topology

	ComputeEdgeLengths(mesh);				// edge lengths
	ComputeElementAreas(mesh);				// element areas
	computeBarryPoints(mesh);				// barry points for each element

	extractInterface(mesh);					// interface nodes and edges

	/* 5. Plot mesh */
	if (is_to_plot)
	{	gnuplot_mesh(mesh, filename+".1");
		gnuplot_interface(mesh, filename+".1");
	}


	return mesh;
}




}
