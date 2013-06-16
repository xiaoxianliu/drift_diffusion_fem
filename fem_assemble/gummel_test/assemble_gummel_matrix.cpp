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

	string filename="assemble_gummel_mat";

	/* 1. Generate mesh */
	MeshData mesh;

	// interface nodes
	vector<vector<double> > interface_nodes;
	vector<double> new_node(2);
	new_node[0] = 0.0;	new_node[1]=1.0;
	interface_nodes.push_back(new_node);
	new_node[0] = 0.0;	new_node[1]=0.0;
	interface_nodes.push_back(new_node);

	// generate mesh
	mesh = generateMesh(filename, interface_nodes);


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
