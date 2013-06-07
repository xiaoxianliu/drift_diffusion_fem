#include <vector>

#include <armadillo>
#include "../triangle/mesh.hpp"
#include "../fem_assemble/fem_assemble.hpp"

#include "exciton_cos_interface.hpp"

#include "misc.hpp"			// include functions related to the definition of state and adjoint equations


/* Main solver for state equation */
arma::vec solveStateEq(const my_mesh::MeshData mesh)
{
using namespace arma;
	/* 1. Assemble linear system (without imposing Dirichlet BC)*/
	/*	M*u = g_vec	*/

	/* 1.1 Coefficient matrix */
	arma::mat M;
	{	arma::vec vec_a = linear_fem::interpolateFunction(mesh, func_a);
		arma::mat A = linear_fem::assembleMatrixA(mesh, vec_a);

		arma::vec vec_c = linear_fem::interpolateFunction(mesh, func_c);
		arma::mat C = linear_fem::assembleMatrixC(mesh, vec_c);

		arma::vec vec_d = linear_fem::interpolateFunction(mesh, func_d);
		arma::mat D = linear_fem::assembleMatrixD(mesh, vec_d);

		M = A + C + D;
	}

	/* 1.2 right-hand side vector */
	arma::vec g_vec = linear_fem::interpolateFunction(mesh, func_g);

	/* 2. Apply Dirichlet boundary condition */
	int num_nodes = mesh.num_nodes;
	for (int i=0; i<num_nodes; i++)
	{
		std::vector<int> Es = mesh.topology0to1[i];			// all neighboring edges to i-th vertex
		for (int j=0; j<Es.size(); j++)
		{	int edge_index = Es[j];
			if ( mesh.edge_markers[edge_index] == 1 || mesh.edge_markers[edge_index] == 3 )
			{	double x = mesh.nodes[i][0];
				double y = mesh.nodes[i][1];
				M(i, arma::span::all) = arma::zeros<arma::mat>(1, num_nodes);
				M(i,i) = 1.0;
				g_vec(i) = uD(x,y);

				break;
			}
		}
	}

	/* 3. Solve for solution */
	arma::vec u;
	u = arma::solve(M, g_vec);

	return u;
}


