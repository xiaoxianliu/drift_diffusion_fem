#include <vector>
#include <armadillo>
#include "../triangle/mesh.hpp"
#include "../fem_assemble/fem_assemble.hpp"
#include "inorganicDD.hpp"

#include "parameters.hpp"

// local function declaration
int applyDirichletBC_n(	const my_mesh::MeshData &mesh, const arma::vec n_D, arma::mat Mn, arma::vec rhs);


/** Main solver for NContinuity equation */
arma::vec solveNContinuityEq(	const my_mesh::MeshData &mesh,
				const arma::vec &input_psi,
				const arma::vec &input_n,
				const arma::vec &input_p,
				const arma::vec &n_D)				// Dirichlet boundary condition
{
using namespace linear_fem;


	// 1. Assemble coefficient matrix
	arma::vec mu_n_vec = interpolateConstant(mesh, parameters::mu_n);
	arma::mat Mn = assembleMatrixGummel(mesh, input_psi, mu_n_vec);

	// 2. Assemble right-hand-side vector
	arma::vec rhs;
	arma::vec g_vec = interpolateFunction(mesh, generation_Func);
	arma::vec recomb_vec = recombSRH_Vec(input_n, input_p);
	rhs = L2project_Vec(mesh, g_vec-recomb_vec);

	// 3. Apply Dirichlet boundary conditions
	applyDirichletBC_n(mesh, n_D, Mn, rhs);

	// 4. Solve for and return solution
	arma::vec output_n;
	output_n = arma::solve(Mn, rhs);

	return output_n;
}


int applyDirichletBC_n(	const my_mesh::MeshData &mesh, const arma::vec n_D, arma::mat Mn, arma::vec rhs)
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
		{	Mn(i, arma::span::all).zeros();
			Mn(i,i) = 1.0;
			rhs(i) = n_D(i);
		}
	}

	return 0;
}
