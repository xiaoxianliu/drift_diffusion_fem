#include <iostream>
#include <cstdlib>

#include <armadillo>

#include "../triangle/mesh.hpp"
#include "../my_fem/my_fem.hpp"
#include "auxillary.hpp"

/************************************************************************************************************************************/
// Electric field: 	dE/dnu = -dpsi/dnu
// (1) On Gamma_anode:	dE/dnu1
// (2) On Gamma_interface:	dE/dnu1, dE/dnu2
// (3) On Gamma_cathod:	dE/dnu2

/************************************************************************************************************************************/
// Normal component of Electron flux:		Fn_dot_nu
// (1) On Gamma_anode:	Fn_dot_nu1

int compute_Boundary1Flux_n(	const my_mesh::MeshData &mesh,
				const arma::vec &n,
				const arma::vec &psi,
				arma::vec &Fn_nu1)
{
	int num_boundary_nodes = mesh.boundary1_nodes.size();

	// 0. Compute relevant quantities
	// Electric field
	arma::vec E(mesh.num_nodes);
	compute_ElectricFieldAmplitude(mesh, psi, E);
	// piecewise mobility for electron
	arma::vec mu_n_elementwise(mesh.num_elements);
	compute_MobilityN_elementwise(mesh, E, mu_n_elementwise);
	
	// 1. Coefficient matrix
	arma::mat M(num_boundary_nodes, num_boundary_nodes);
	M.zeros();
	for (int i=0; i<mesh.boundary1_edges.size(); i++)			// loop through boundary edges
	{	int edge_index = mesh.boundary1_edges[i];
		double length = mesh.edge_lengths [edge_index];
		M(i,i) += 0.5 * length;
		M(i+1, i+1) += 0.5 * length;
	}

	// 2. right-hand-side vector
	arma::vec rhs(num_boundary_nodes);
	arma::mat M_gummel = my_fem::assembleMatrixGummel(mesh, -psi, mu_n_elementwise);
	for (int i=0; i<num_boundary_nodes; i++)
	{	int node_index = mesh.boundary1_nodes[i];			// global index of current node

		arma::vec M_gummel_times_n_vec = M_gummel * n;
		rhs(i) = - M_gummel_times_n_vec(node_index);
	}

	// 3. solve for normal electron flux "Fn" on boundary 1
	Fn_nu1 = arma::solve(M, rhs);

return 0;
}
// (2) On Gamma_interface:	Fn_dot_nu1, Fn_dot_nu2
// (3) On Gamma_cathode:	Fn_dot_nu2


/************************************************************************************************************************************/
// Normal component of hole flux:		Fp_dot_nu
// (1) On Gamma_anode:	Fp_dot_nu1

int compute_Boundary1Flux_p(	const my_mesh::MeshData &mesh,
				const arma::vec &p,
				const arma::vec &psi,
				arma::vec &Fp_nu1)
{
	int num_boundary_nodes = mesh.boundary1_nodes.size();

	// 0. Compute relevant quantities
	// Electric field
	arma::vec E(mesh.num_nodes);
	compute_ElectricFieldAmplitude(mesh, psi, E);
	// piecewise mobility for electron
	arma::vec mu_p_elementwise(mesh.num_elements);
	compute_MobilityP_elementwise(mesh, E, mu_p_elementwise);
	
	// 1. Coefficient matrix
	arma::mat M(num_boundary_nodes, num_boundary_nodes);
	M.zeros();
	for (int i=0; i<mesh.boundary1_edges.size(); i++)			// loop through boundary edges
	{	int edge_index = mesh.boundary1_edges[i];
		double length = mesh.edge_lengths [edge_index];
		M(i,i) += 0.5 * length;
		M(i+1, i+1) += 0.5 * length;
	}

	// 2. right-hand-side vector
	arma::vec rhs(num_boundary_nodes);
	arma::mat M_gummel = my_fem::assembleMatrixGummel(mesh, psi, mu_p_elementwise);
	for (int i=0; i<num_boundary_nodes; i++)
	{	int node_index = mesh.boundary1_nodes[i];			// global index of current node

		arma::vec M_gummel_times_p_vec = M_gummel * p;
		rhs(i) = - M_gummel_times_p_vec(node_index);
	}

	// 3. solve for normal electron flux "Fn" on boundary 1
	Fp_nu1 = arma::solve(M, rhs);

return 0;
}
// (2) On Gamma_interface:	Fp_dot_nu1, Fp_dot_nu2
// (3) On Gamma_cathode:	Fp_dot_nu2

/************************************************************************************************************************************/
// Normal component of exciton flux:		Fx_dot_nu
// (1) On Gamma_anode:	Fx_dot_nu1
// (2) On Gamma_interface:	Fx_dot_nu1, Fx_dot_nu2
// (3) On Gamma_cathode:	Fx_dot_nu2
