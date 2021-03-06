#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include <armadillo>

#include "../../triangle/mesh.hpp"
#include "../../my_fem/my_fem.hpp"
#include "parameters.hpp"
#include "auxillary.hpp"

// Solver for the continuity equation of holes 


/***************************************************************************************************************************/
// Declaration of local functions
int apply_DirichletBC_p(const my_mesh::MeshData &mesh,			// apply Dirichlet boundary conditions
			arma::mat &M,
			arma::vec &rhs);



/***************************************************************************************************************************/
// Main function to solve hole continuity equation in Gummel's iteration
int solve_ContinuityEq_p(const my_mesh::MeshData &mesh, 
			const arma::vec &input_psi,
			const arma::vec &input_n,
			const arma::vec &input_p,
			const arma::vec &input_u,
			arma::vec &output_p)
{
	// 1. compute related functions/vectors
	// 1.1 electric field amplitude
	arma::vec E(mesh.num_nodes);
	compute_ElectricFieldAmplitude(mesh, input_psi, E);
	// 1.2 element-wise mobility
	arma::vec mu_p_elementwise(mesh.num_elements);
	compute_MobilityP_elementwise (mesh, E, mu_p_elementwise);
	// 1.3 recombination rate
	arma::vec recomb_interface;
	compute_RecombinationRate(mesh, E, recomb_interface);
	// 1.4 Exciton dissociation rate
	arma::vec k_diss_interface;
	compute_ExcitonDissociationRate(mesh, input_psi, input_n, input_p, input_u, k_diss_interface);

	// 2. Assemble coefficient matrix
	arma::mat M;
	{	// 2.1 assemble Gummel matrix
		arma::mat M1 = my_fem::assembleMatrixGummel(mesh, input_psi, mu_p_elementwise);
		// 2.2 "matrixD" from the interface integral of "p"
		arma::mat M2 = my_fem::assembleMatrixD(mesh, recomb_interface % input_n);

		M = M1 + M2;
	}

	// 3. Assemble right-hand-side vector
	arma::vec rhs;
	{	arma::mat M_rhs = my_fem::assembleMatrixD(mesh, k_diss_interface);
		rhs = M_rhs * input_u;
	}

	// 4. Apply Dirichlet boundary conditions
	apply_DirichletBC_p(mesh, M, rhs);

	// 5. Solve for solution "output_p"
	output_p = arma::solve(M, rhs);


return 0;
}





/***************************************************************************************************************************/
// Definition of local functions
// Dirichlet boundary condition for p-continuity equation
int apply_DirichletBC_p(const my_mesh::MeshData &mesh,
			arma::mat &M,
			arma::vec &rhs)
{
	for (int i=0; i<mesh.num_nodes; i++)
	{
		// 1. Determine if node "i" is on anode or cathode
		std::vector<int> neigh_edges = mesh.topology0to1[i];		// identify neighboring edges
		bool is_on_anode = false, is_on_cathode = false;	// anode is left boundary "1"; cathode is right boundary "3"

		for (int j=0; j<neigh_edges.size(); j++)		
		{	int edge_index = neigh_edges[j];
			if (mesh.edge_markers[edge_index] == 1)
			{	is_on_anode=true;	break;	}
			else if (mesh.edge_markers[edge_index] == 3)
			{	is_on_cathode = true;	break;	}
		}

		// 2. If is on anode or cathod, set the corresponding boundary values
		if (is_on_anode)
		{	M(i, arma::span::all) = arma::zeros<arma::mat>(1, mesh.num_nodes);
			M(i,i) = 1.0;
			rhs(i) = parameters::p_A;
		}
		else if (is_on_cathode)
		{	M(i, arma::span::all) = arma::zeros<arma::mat>(1, mesh.num_nodes);
			M(i,i) = 1.0;
			rhs(i) = parameters::p_C;
		}
	}

return 0;
}
