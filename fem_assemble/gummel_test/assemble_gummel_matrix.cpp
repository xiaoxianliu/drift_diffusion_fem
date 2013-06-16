#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include "../fem_assemble.hpp"
#include "../../triangle/mesh.hpp"



int main(int argc, char* argv[])
{
using namespace std;
using namespace arma;
using namespace my_mesh;
using namespace linear_fem;


	/* 1. Generate mesh */
	MeshData mesh;

	// initialize "mesh"
	vector<vector<double> > nodes;
	vector<int> node_markers;
	vector<vector<int> > segments;
	vector<int> segment_markers;
	vector<vector<double> > regions;
	vector<int> region_markers;


	vector<double> new_node(2);
	new_node[0] = 0.0;	new_node[1]=1.0;
	nodes.push_back(new_node);
	new_node[0] = 0.0;	new_node[1]=0.0;
	nodes.push_back(new_node);


	newMesh(nodes, node_markers,		// 4 corner vertices are added to "nodes"
		segments, segment_markers,
		regions, region_markers);

	string filename = "test_gummel";
	writePolyfile(filename, nodes, node_markers, segments, segment_markers, regions, region_markers);

	// run "triangle" to generate a coarse mesh
	string cmd = "/home/xiaoxian/bin/triangle/triangle -qzpeAa1 " + filename + ".poly";
	system(cmd.c_str());

	string node_file_name= filename + ".1.node";
	string edge_file_name= filename + ".1.edge";
	string ele_file_name = filename + ".1.ele";

	ReadNodes(mesh, node_file_name);
	ReadEdges(mesh, edge_file_name);
	ReadElements(mesh, ele_file_name);

	ComputeTopology(mesh);

	ComputeEdgeLengths(mesh);
	ComputeElementAreas(mesh);

	gnuplot_mesh(mesh, filename);


	/* 2. Generate Gummel matrix and write it to output */
	/*	When potential "psi" is constant, gradient of "psi" is zero vector. 
	/*	Hence, if A is the matrix of  integral[(grad(u), grad(v)) ], one should have
			M_gummel == A
	*/
	mat A, M_gummel;

	vec psi(mesh.num_nodes), mu(mesh.num_elements);
	psi.ones();	mu.ones();
	M_gummel = assembleMatrixGummel(mesh, psi, mu);

	vec a(mesh.num_nodes);
	a.ones();
	A = assembleMatrixA(mesh, a);

	cout << M_gummel << "\n";
	cout << A << "\n";
	cout << M_gummel - A << "\n";

	
	return 0;
}
