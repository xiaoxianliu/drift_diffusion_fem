#include <cstdlib>
#include <string>
#include "mesh.hpp"

using namespace std;

int TriangleMesh(MeshData &mesh, string& filename){
	/* Generate mesh by "triangle" from an initial input file*/
	string cmd = "/home/xiaoxian/bin/triangle/triangle -qzpeAa0.0005 "; 
	if (filename.length()==0)
		filename = "rectangle0";
	cmd += (filename + ".poly");
	system(cmd.c_str());

	/* Read in mesh information and compute mesh quantities */
	string node_file_name= filename + ".1.node";
	string edge_file_name= filename + ".1.edge";
	string ele_file_name = filename + ".1.ele";

	ReadNodes(mesh, node_file_name);
	ReadEdges(mesh, edge_file_name);
	ReadElements(mesh, ele_file_name);

	/* Compute mesh topology */
	ComputeTopology(mesh);

	/* Compute geometric quantities */
	ComputeEdgeLengths(mesh);
	ComputeElementAreas(mesh);

	/* Plot with GNU-plot */
	WriteGNUplot(mesh, filename);

	return 0;
}
