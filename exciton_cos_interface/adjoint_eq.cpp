#include <vector>

#include <armadillo>

#include "../triangle/mesh.hpp"
#include "../fem_assemble/fem_assemble.hpp"
#include "exciton_cos_interface.hpp"

double func_a_adjoint(double x, double y)
{	return 1.0;
}
double func_c_adjoint(double x, double y)
{	return 1.0;
}
double func_d_adjoint(double x, double y)
{	return 1.0;
}
double xiD(double x, double y)				// xi_D is Dirichlet boundary condition
{	return 0.0;
}



/* Main solver for state equation */
arma::vec solveAdjointEq(const my_mesh::MeshData &mesh)
{

	int num_nodes = mesh.num_nodes;

	/* 1. Assemble coefficient matrix */
	/*	M * xi = rhs	*/

	/* 1.1 Coefficient matrix */
	arma::mat M;
	{
	arma::vec vec_a = linear_fem::interpolateFunction(mesh, func_a_adjoint) ;
	arma::mat A = linear_fem::assembleMatrixA(mesh, vec_a);

	arma::vec vec_c = linear_fem::interpolateFunction(mesh, func_c_adjoint);
	arma::mat C = linear_fem::assembleMatrixC(mesh, vec_c);

	arma::vec vec_d = linear_fem::interpolateFunction(mesh, func_d_adjoint);
	arma::mat D = linear_fem::assembleMatrixD(mesh, vec_d);

	M = A + C + D;
	}
	/* 1.2 right-hand side vector */
	arma::vec rhs_vec;
	{
	arma::vec vec_ones = arma::ones<arma::vec>(num_nodes);
	arma::mat D = linear_fem::assembleMatrixD(mesh, vec_ones);
	rhs_vec = D*vec_ones;
	}

	/* 2. apply Dirichlet boundary condition */
	for (int i=0; i<num_nodes; i++)
	{
		std::vector<int> Es = mesh.topology0to1[i];			// all neighboring edges to i-th vertex
		for (int j=0; j<Es.size(); j++)
		{	int edge_index = Es[j];
			if ( mesh.edge_markers[edge_index] == 1 || mesh.edge_markers[edge_index] == 3 )
			{	double x = mesh.nodes[i][0];
				double y = mesh.nodes[i][1];

				M(i, arma::span::all) = arma::zeros<arma::mat>(1, num_nodes);	// modify matrix on Dirichlet node
				M(i,i) = 1.0;
				rhs_vec(i) = xiD(x,y);						// modify rhs vector on Dirichlet node

				break;
			}
		}
	}


	/* 3. Solve for solution */
	arma::vec xi;
	xi = arma::solve(M, rhs_vec);


	return xi;
}
