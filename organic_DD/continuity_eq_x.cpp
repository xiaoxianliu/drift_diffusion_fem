#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <armadillo>

#include "../triangle/mesh.hpp"
#include "../my_fem/my_fem.hpp"

#include "parameters.hpp"
#include "auxillary.hpp"

// Solving for exciton density "u" for given {input_psi, input_n, input_p, input_u} in Gummel's iteration


/******************************************************************************************************************/
// Declarations of local functions

int apply_DirichletBC_x(const my_mesh::MeshData &mesh,
			arma::mat &M,
			arma::vec &rhs);





/******************************************************************************************************************/
// Main function

int solve_ContinuityEq_x(	const my_mesh::MeshData &mesh,
				const arma::vec &input_psi,
				const arma::vec &input_n,
				const arma::vec &input_p,
				const arma::vec &input_u,
				arma::vec &output_u)
{
	// 1. Compute related functions/vectors
	// 1.1 vector of exciton mobility (constant vector);
	arma::vec mu_x_vec(mesh.num_nodes);
	mu_x_vec.ones();
	mu_x_vec = mu_x_vec * parameters::mu_x;
	// 1.2 vector of exciton decay rate (1/exciton_lifetime)
	arma::vec decay_rate_x(mesh.num_nodes);
	decay_rate_x.ones();
	decay_rate_x = decay_rate_x * (1.0 / parameters::tau_x);
	// 1.3 vector of interface dissociation rate of exciton into free carriers
	arma::vec k_diss_interface;
	compute_ExcitonDissociationRate (mesh, input_psi, input_n, input_p, input_u, k_diss_interface);
	// 1.4 vector of photo generation rate of exciton
	arma::vec Q_vec(mesh.num_nodes);
	compute_PhotoGenerationVec(mesh, Q_vec);
	// 1.5 vector of electric field amplitude
	arma::vec E;
	compute_ElectricFieldAmplitude(mesh, input_psi, E);
	// 1.6 vector of interface recombination rate of free carriers into excitons
	arma::vec recomb_interface;
	compute_RecombinationRate(mesh, E, recomb_interface);


	// 2. Assemble coefficient matrix
	arma::mat M;
	{
		// 2.1 matrix from diffusion operator
		arma::mat M1 = my_fem::assembleMatrixA(mesh, mu_x_vec);
		// 2.2 matrix from exciton decay
		arma::mat M2 = my_fem::assembleMatrixC(mesh, decay_rate_x);
		// 2.3 matrix from interface integral of exciton dissociation
		arma::mat M3 = my_fem::assembleMatrixD(mesh, k_diss_interface);

		M = M1 + M2 + M3;
	}

	// 3. Assemble right-hand-side vector
	arma::vec rhs;
	{
		// 3.1 vector from volume integral of photo generation
		arma::mat M_rhs1 = my_fem::assembleMatrixC(mesh, arma::ones<arma::vec>(mesh.num_nodes));
		arma::vec rhs1 = M_rhs1 * Q_vec;
		// 3.2 vector from interface integral of recombination from free carriers into excitons
		arma::mat M_rhs2 = my_fem::assembleMatrixD(mesh, recomb_interface);
		arma::vec rhs2 = M_rhs2 * (input_n % input_p);

		rhs = rhs1 + rhs2;
	}

	// 4. Apply Dirichlet boundary conditions
	apply_DirichletBC_x(mesh, M, rhs);

	// 5. Solve for new exciton density
	output_u = arma::solve(M, rhs);

return 0;
}










/******************************************************************************************************************/
// Definitions of local functions



int apply_DirichletBC_x(const my_mesh::MeshData &mesh,
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
			rhs(i) = parameters::u_A;
		}
		else if (is_on_cathode)
		{	M(i, arma::span::all) = arma::zeros<arma::mat>(1, mesh.num_nodes);
			M(i,i) = 1.0;
			rhs(i) = parameters::u_C;
		}
	}

return 0;
}

