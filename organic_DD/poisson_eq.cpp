#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <armadillo>

#include "../triangle/mesh.hpp"
#include "../my_fem/my_fem.hpp"

#include "parameters.hpp"
#include "auxillary.hpp"


// Solving for solution of nonlinear Poisson's equation


/**********************************************************************************************************/
// Declarations of local functions

// Solver for linearized Poisson's equation
int solve_LinearizedPoissonEq(	const my_mesh::MeshData &mesh,
				const arma::vec &input_psi,
				const arma::vec &input_n,
				const arma::vec &input_p,
				const arma::vec &input_u,
				arma::vec &output_psi,
				double applied_psi);







/**********************************************************************************************************/
// Main function
int solve_NonlinearPoissonEq(	const my_mesh::MeshData &mesh,
				const arma::vec &input_psi,
				const arma::vec &input_n,
				const arma::vec &input_p,
				const arma::vec &input_u,
				arma::vec &output_psi,
				double applied_psi)
{

	// 1. Set up Newton's iteration of solving nonlinear Poisson's equation
	//	Note: when max_iter = 1, effectively one recovers the Gummel's iteration
	int max_iter=1;
	double newton_err=1.0, newton_tol=1e-3;
	arma::vec prev_psi = input_psi;				// previous "psi" for each iteration
	arma::vec new_psi;					// new "psi" solved for after each iteration

	// 2. Newton's iteration
	for (int iter=0; iter < max_iter; iter++)
	{
//		std::cout << "\tIn Newton's solver for Poisson's equation, iter " << iter << "\n";
		// Solve linearized Poisson's equation
		solve_LinearizedPoissonEq(mesh, prev_psi, input_n, input_p, input_u, new_psi, applied_psi);

		// Compute Newton's error
		newton_err = arma::norm(prev_psi-new_psi, "inf");
		if (newton_err < newton_tol)
			break;
		prev_psi = new_psi;
	}

	output_psi = new_psi;

return 0;
}









/**********************************************************************************************************/
// Definition of local functions


// Dirichlet boundary condition for linearized Poisson's equation
int apply_DirichletBC_LinearizedPoissonEq(	const my_mesh::MeshData &mesh,
						arma::mat &M,
						arma::vec &rhs,
						double applied_psi)		// applied potential
{
	for (int i=0; i< mesh.num_nodes; i++)
	{
		// 1. Determine if node "i" is on the Dirichlet boundary: "1" --> anode, "2" --> cathode
		bool is_on_anode=false, is_on_cathode=false;
		std::vector<int> neigh_edges = mesh.topology0to1[i];
		for (int j=0; j<neigh_edges.size(); j++)
		{	int edge_index = neigh_edges[j];

			if (mesh.edge_markers[edge_index]==1)
			{	is_on_anode = true;	break;
			}
			else if (mesh.edge_markers[edge_index]==3)
			{	is_on_cathode = true;	break;
			}
		}

		// 2. Modify the entries in "M" and "rhs" for the Dirichlet nodes
		if (is_on_anode)
		{	M(i, arma::span::all) = arma::zeros<arma::mat>(1, mesh.num_nodes);
			M(i,i)=1.0;
			rhs(i) = parameters::psi_A + applied_psi;		// applied_psi is always assumed on anode!!!!
		}
		else if (is_on_cathode)
		{	M(i, arma::span::all) = arma::zeros<arma::mat>(1, mesh.num_nodes);
			M(i,i)=1.0;
			rhs(i) = parameters::psi_C;
		}
	}

return 0;
}



















// Solve linearized Poisson's equation for each Newton's iteration
int solve_LinearizedPoissonEq(	const my_mesh::MeshData &mesh,
				const arma::vec &input_psi,
				const arma::vec &input_n,
				const arma::vec &input_p,
				const arma::vec &input_u,
				arma::vec &output_psi,
				double applied_psi)
{
using parameters::lambda_squared;
using parameters::epsilon_rel;

	// 1. Assemble coefficient matrix
	arma::mat M(mesh.num_nodes, mesh.num_nodes);
	{
		// 2.1 coefficient of laplacian operator
		arma::vec a1_vec = arma::ones<arma::vec>(mesh.num_nodes) * lambda_squared * epsilon_rel;
		arma::mat M1 = my_fem::assembleMatrixA(mesh, a1_vec);

		arma::mat M2 = my_fem::assembleMatrixC(mesh, input_p + input_n);

		M = M1 + M2;
	}

	// 2. Assemble right-hand side vector
	arma::vec rhs;
	{
		arma::vec rhs1 = my_fem::L2project_Vec(mesh, input_psi % (input_p + input_n));
		arma::vec rhs2 = my_fem::L2project_Vec(mesh, input_p - input_n);
		rhs = rhs1 + rhs2;
	}

	// 3. Apply Dirichlet boundary conditions
	apply_DirichletBC_LinearizedPoissonEq(mesh, M, rhs, applied_psi);

	// 4. Solve for solution
	output_psi = arma::solve(M, rhs);
return 0;
}


