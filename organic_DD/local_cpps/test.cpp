#include <iostream>
#include <vector>

#include <armadillo>


#include "../../triangle/mesh.hpp"
#include "../../my_fem/my_fem.hpp"
#include "../../gnuplot/my_gnuplot.hpp"

#include "../parameters.hpp"

/*******************************************************************************************************************/
// Test function that solves "-lambda_squared * grad( epsilon_r*grad(psi) ) = p - n" for given "p" and "n"
// purpose: test the consistency of solution "psi, n, p" from Gummel's iteration


// Local function: Dirichlet boundary conditions
int apply_DirichletBC_psi(	const my_mesh::MeshData &mesh, 
				arma::mat &M,
				arma::vec &rhs,
				double applied_psi)
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






// Main function: 
int test_Compute_psi(	const my_mesh::MeshData &mesh,
			const arma::vec &n,
			const arma::vec &p,
			arma::vec &psi,
			double applied_psi)
{
using parameters::lambda_squared;
using parameters::epsilon_rel;

	// 1. Assemble coefficient matrix
	arma::vec a1_vec = arma::ones<arma::vec>(mesh.num_nodes) * lambda_squared * epsilon_rel;
	arma::mat M = my_fem::assembleMatrixA (mesh, a1_vec);

	// 2. Assemble right-hand-side vector
	arma::vec rhs = my_fem::L2project_Vec(mesh, p-n);

	// 3. Apply Dirichlet boundary condition
	apply_DirichletBC_psi(mesh, M, rhs, applied_psi);

	// 4. solve for solution
	psi = arma::solve(M, rhs);
return 0;
}




/***********************************************************************************************************************/
