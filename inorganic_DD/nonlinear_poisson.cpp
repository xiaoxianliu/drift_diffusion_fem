#include <vector>
#include <armadillo>
#include "../triangle/mesh.hpp"
#include "../fem_assemble/fem_assemble.hpp"

#include "parameters.hpp"
//#include "inorganicDD.hpp"





/*************       Dirichlet boundary condition when linearized Poisson's equation is solved **************/
int applyDirichletBC_psi(	const my_mesh::MeshData &mesh,
				const arma::vec &psi_D,
				arma::mat &M,
				arma::vec &rhs)
{
	for (int i=0; i<mesh.num_nodes; i++)
	{
		// 1. Determine if node "i" is on Dirichlet boundaries
		bool is_on_dirichlet_boundary = false;
		std::vector<int> neigh_edges = mesh.topology0to1[i];
		for (int j=0; j<neigh_edges.size(); j++)
		{	int edge = neigh_edges[j];
			if (mesh.edge_markers[edge] == 1 || mesh.edge_markers[edge] == 3)
			{	is_on_dirichlet_boundary = true;	break;	}
		}

		// 2. If it's on the Dirichlet boundaries, modify coefficient matrix and right-hand-side vector
		if (is_on_dirichlet_boundary)
		{	M(i, arma::span::all).zeros();
			M(i,i) = 1.0;
			rhs(i) = psi_D(i);
		}
	}
	return 0;
}



/**************     Linearized poisson equation for each Newton's iteration ******************/
arma::vec solveLinearizedPoissonEq(	const my_mesh::MeshData &mesh,
					const arma::vec &input_psi,
					const arma::vec &input_n,
					const arma::vec &input_p,
					const arma::vec &psi_D,
					const arma::vec &C	)			// "C" is arma::vec form of doping density
{
using namespace linear_fem;

	// 1. Assemble coefficient matrix
	arma::mat M;
	{	arma::vec a_vec = interpolateConstant(mesh, parameters::lambda_squared);
		arma::mat A = assembleMatrixA(mesh, a_vec);			// int [ a < grad(u), grad(v) > ]
		arma::vec c_vec = input_n + input_p;
		arma::mat C = assembleMatrixC(mesh, c_vec);			// int [ cuv ]
		M = A + C;
	}

	// 2. Assemble right-hand-side matrix
	arma::vec dummy_vec = (input_n + input_p) % input_psi + C + input_p - input_n;	// element-wise formula for rhs function
	arma::vec rhs = L2project_Vec(mesh, dummy_vec);

	// 3. Apply Dirichlet boundary condition
	applyDirichletBC_psi(mesh, psi_D, M, rhs);

	// 4. Solve linear equations
	arma::vec output_psi = arma::solve(M, rhs);

	return output_psi;
}











/****************    Main function:  Newton's method solving nonlinear Poisson's equation      ******************/
arma::vec solveNonlinearPoissonEq(	const my_mesh::MeshData &mesh,
					const arma::vec &input_psi,
					const arma::vec &input_n,
					const arma::vec &input_p,
					const arma::vec &psi_D,
					const arma::vec &C	)			// "C" is arma::vec form of doping density
{
	arma::vec prev_psi = input_psi;			// initial "psi" to start Newton's iteration
	arma::vec new_psi;

	double newton_err=1.0, newton_tol = 1e-3;
	int max_iter = 10;
	for (int iter; iter < max_iter; iter++)
	{
		new_psi = solveLinearizedPoissonEq(mesh, input_psi, input_n, input_p, psi_D, C);
		newton_err = arma::norm( arma::abs( new_psi - prev_psi ), "inf");
		if (newton_err < newton_tol)	break;
		else
		{	prev_psi = new_psi;	}
	}

	return new_psi;
}
