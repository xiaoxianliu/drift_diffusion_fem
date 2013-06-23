#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

#include "exciton_smart_geometry.hpp"
#include "../my_fem/my_fem.hpp"

int main()
{
using namespace std;

	std::string filename = "smart";

	/* 1. Generate mesh */

	vector<vector<double> > interface_nodes;
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

	double max_area = 0.001;
	my_mesh::MeshData mesh = my_mesh::generateMesh(filename, interface_nodes, max_area);


	/* 2. Form coefficient matrix and RHS vector */
	/* 2.1 Coefficient matrix */
	arma::vec a = interpolate_func(mesh, diffusion_coefficient);
	arma::vec c = interpolate_func(mesh, decay_rate);
	arma::vec d = interpolate_func(mesh, interface_reaction_rate);

	arma::mat A = my_fem::assembleMatrixA(mesh, a);
	arma::mat C = my_fem::assembleMatrixC(mesh, c);
	arma::mat D = my_fem::assembleMatrixD(mesh, d);

	arma::mat M = A + C + D;

	/* 2.2 RHS vector */
	arma::vec g = project_func(mesh, generation_rate);

	/* 2.3 Apply Dirichlet boundary conditions to both coefficient matrix and rhs vector */
	applyDirichletBC(mesh, M, g);

	/* 3 Solve for and plot solution */
	/* 3.1 solve for solution */
	arma::vec u;
	u = arma::solve(M,g);

	/* 3.2 plot solution */
	plotSolution(mesh, u, filename);


	return 0;
}



