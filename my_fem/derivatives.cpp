#include <vector>
#include <string>
#include <fstream>
#include <armadillo>
#include "../triangle/mesh.hpp"
#include "my_fem.hpp"


namespace my_fem
{

int computeGradient(	const my_mesh::MeshData &mesh,
			const arma::vec &psi,
			arma::vec &psi_grad_x,
			arma::vec &psi_grad_y)
{

	// 0. Check vector dimentions
	int num_nodes = mesh.num_nodes;
	if ( psi_grad_x.n_rows != num_nodes || psi_grad_x.n_cols != 1)
	{	psi_grad_x.resize(num_nodes, 1);	psi_grad_x.zeros();	}
	if ( psi_grad_y.n_rows != num_nodes || psi_grad_y.n_cols != 1)
	{	psi_grad_y.resize(num_nodes, 1);	psi_grad_y.zeros();	}


	// 1. Assemble coefficient matrix for Galerkin method
	arma::mat M = my_fem::assembleMatrixC (mesh, arma::ones<arma::vec>(num_nodes)) ;

	// 2. Assemble rhs vector of functionals for both psi_grad_x and psi_grad_y, respectively
	arma::vec rhs_x(num_nodes);	rhs_x.zeros();
	arma::vec rhs_y(num_nodes);	rhs_y.zeros();

	for (int t=0; t<mesh.num_elements; t++)
	{
		int v0 = mesh.topology2to0[t][0];				// indices of 3 nodes of triangle "t"
		int v1 = mesh.topology2to0[t][1];
		int v2 = mesh.topology2to0[t][2];

		arma::mat Jacobian(2,2);					// Jacobian matrix of dx/dxi ("xi" is coordinates of
		Jacobian(0,0) = mesh.nodes[v1][0] - mesh.nodes[v0][0];		//   reference triangle
		Jacobian(1,0) = mesh.nodes[v1][1] - mesh.nodes[v0][1];
		Jacobian(0,1) = mesh.nodes[v2][0] - mesh.nodes[v0][0];
		Jacobian(1,1) = mesh.nodes[v2][1] - mesh.nodes[v0][1];

		arma::vec dpsi_dxi(2);						// dpsi / dxi on the reference triangle
		dpsi_dxi(0) = psi(v1) - psi(v0);
		dpsi_dxi(1) = psi(v2) - psi(v0);

		arma::vec psi_gradient(2);
		psi_gradient = Jacobian.i().t() * dpsi_dxi;		// dpsi / dx = Jacobian^(-T) * dpsi/dxi


		double area_t = mesh.ele_areas[t];
		rhs_x(v0) += area_t * psi_gradient(0) / 3.0;
		rhs_x(v1) += area_t * psi_gradient(0) / 3.0;
		rhs_x(v2) += area_t * psi_gradient(0) / 3.0;

		rhs_y(v0) += area_t * psi_gradient(1) / 3.0;
		rhs_y(v1) += area_t * psi_gradient(1) / 3.0;
		rhs_y(v2) += area_t * psi_gradient(1) / 3.0;

	}

	// 3. Solve for psi_grad_x and psi_grad_y
	psi_grad_x = solve(M, rhs_x);
	psi_grad_y = solve(M, rhs_y);


	return 0;
}



}
