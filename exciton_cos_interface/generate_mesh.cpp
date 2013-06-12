#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../triangle/mesh.hpp"

#define PI 3.1415926535897932384626433832795

my_mesh::MeshData generateMesh_CosInterface(const std::string &filename, double y_control)
{
using namespace std;
using namespace my_mesh;

	MeshData mesh;


	/* 0. Determine if "y_control" is valid */
	double y_control_max = 0.9;
	if (std::fabs(y_control) > y_control_max)
	{	std::cout << "abs(y_control) has to be smaller than " << y_control_max << "\n";
		exit(1);
	}


	/* 1. Generate sample points on interface */
/* all sample points with coordinates (x,y) satisfy the formula:
x = y_control * (1 - cos(2*pi*y)); */

	/* 1.1 Inteface node coordinates */
	int num_interface_nodes = 201;
	vector< vector<double> > interface_nodes;

	for (int i=0; i<num_interface_nodes; i++)
	{	vector<double> node(2);
		double y;
		if (i<num_interface_nodes-1)
			y = 1.0 - i/static_cast<double>(num_interface_nodes-1);
		else
			y = 0.0;
		double x = y_control * 0.5*(1 - std::cos(2*PI*y));

		node[0] = x;
		node[1] = y;

		interface_nodes.push_back(node);
	}


	/* 1.2 Form other "newMesh" input */
	vector<int> node_markers;
	vector<vector<int> > segments;
	vector<int> segment_markers;
	vector<vector<double> > regions;
	vector<int> region_markers;

	newMesh(interface_nodes, node_markers,		// 4 corner vertices are added to "interface_nodes" inside "newMesh(...)"
		segments, segment_markers,
		regions, region_markers);

	vector<vector<double> > nodes = interface_nodes;



	/* 2. Write to poly file */
	writePolyfile(	filename,
			nodes, node_markers,
			segments, segment_markers,
			regions, region_markers);




	/* 3. Generate mesh by "triangle" */
	string cmd = "/home/xiaoxian/bin/triangle/triangle -qzpeAa0.001 " + filename + ".poly";
	system(cmd.c_str());




	/* 4. read in mesh information (mesh has been generated in last part) */


	string node_file_name= filename + ".1.node";
	string edge_file_name= filename + ".1.edge";
	string ele_file_name = filename + ".1.ele";

	ReadNodes(mesh, node_file_name);
	ReadEdges(mesh, edge_file_name);
	ReadElements(mesh, ele_file_name);




	/* 5. Compute topology(connectivity) and geometric quantities (areas, lengths) */
	ComputeTopology(mesh);

	/* 6 Extract interface */
	extractInterface(mesh, mesh.interface_edges, mesh.interface_nodes);

	/* 7 Compute geometric quantities: areas, lengths, interface curvature */
	ComputeEdgeLengths(mesh);
	ComputeElementAreas(mesh);
	ComputeInterfaceCurvatures(mesh);


	/* 8. Plot mesh and interface */
	gnuplot_mesh(mesh, filename);
	gnuplot_interface(mesh, filename);

	return mesh;	
}








