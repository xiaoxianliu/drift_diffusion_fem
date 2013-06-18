#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include "../../triangle/mesh.hpp"
#include "../fem_assemble.hpp"

#include "main.hpp"

using namespace std;
using namespace arma;

using namespace my_mesh;
using namespace linear_fem;

int main(int argc, char* argv[])
{
	string filename = "gummel_test";
	// 1. Define mesh
	vector< vector<double> > interface_nodes;
	vector<double> new_node(2);
	new_node[0]=0;	new_node[1]=1;	interface_nodes.push_back(new_node);
	new_node[0]=0;	new_node[1]=0;	interface_nodes.push_back(new_node);

	double max_area=1;
	MeshData mesh = generateMesh(filename, interface_nodes, max_area);


	// 2. Assemble matrix
	// 2.1 coefficient matrix
	vec mu = ones<vec>(mesh.num_nodes);
	vec psi = interpolateFunction(mesh, psi_func);
	mat M = assembleMatrixGummel(mesh, psi, mu);

	// 2.2 right-hand-side vector
	vec rhs = L2project_Func(mesh, f_func);

	// 2.3 Apply Dirichlet boundary conditions
	applyDirichletBC(mesh, M, rhs);

	cout << M << "\n";
	cout << rhs << "\n";

	// 2.4 Solve for solution
	arma::vec u = solve(M, rhs);

	// 2.5 plot solution
	plotSolution (mesh, u, filename);
	
	return 0;
}
