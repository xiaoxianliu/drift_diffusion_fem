#include <cstdlib>
#include <vector>
#include <string>
#include "../triangle/mesh.hpp"






my_mesh::MeshData generateMesh(std::string filename)
{
using namespace my_mesh;
using namespace std;

	/* 1. Form new mesh input */
	vector<vector<double> > interface_nodes;
	vector<int> node_markers;
	vector<vector<int> > segments;
	vector<int> segment_markers;
	vector<vector<double> > regions;
	vector<int> region_markers;

	/* 1.1 Add interface node to the vector of interface nodes */
	vector<double> new_node(2);
	new_node[0] = -0.2;	new_node[1]=1.0;
	interface_nodes.push_back(new_node);
	new_node[0] = -0.2;	new_node[1]=2/3.0;
	interface_nodes.push_back(new_node);
	new_node[0] = 0.3;	new_node[1]=2/3.0;
	interface_nodes.push_back(new_node);
	new_node[0] = 0.3;	new_node[1]=1/3.0;
	interface_nodes.push_back(new_node);
	new_node[0] = -0.2;	new_node[1]=1/3.0;
	interface_nodes.push_back(new_node);
	new_node[0] = -0.2;	new_node[1]=0.0;
	interface_nodes.push_back(new_node);


	newMesh(interface_nodes, node_markers,		// 4 corner vertices are added to "interface_nodes" inside "newMesh(...)"
		segments, segment_markers,
		regions, region_markers);

	vector<vector<double> > nodes = interface_nodes;

	/* 2. Write to poly file */
//	string filename="test_fem";
	writePolyfile(	filename,
			nodes, node_markers,
			segments, segment_markers,
			regions, region_markers);


	/* 3. Generate mesh by "triangle" */
	string cmd = "/home/xiaoxian/bin/triangle/triangle -qzpeAa0.005 " + filename + ".poly";
	system(cmd.c_str());

	/* 4. read in mesh information */
	MeshData mesh;

	string node_file_name= filename + ".1.node";
	string edge_file_name= filename + ".1.edge";
	string ele_file_name = filename + ".1.ele";

	ReadNodes(mesh, node_file_name);
	ReadEdges(mesh, edge_file_name);
	ReadElements(mesh, ele_file_name);

	/* 5. Compute topology(connectivity) and geometric quantities (areas, lengths) */
	ComputeTopology(mesh);

	ComputeEdgeLengths(mesh);
	ComputeElementAreas(mesh);

	/* 6. Plot mesh and interface */
	gnuplot_mesh(mesh, filename);
	gnuplot_interface(mesh, filename);

	return mesh;	
}


