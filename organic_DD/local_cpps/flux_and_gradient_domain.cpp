#include <iostream>
#include <vector>

#include <armadillo>

#include "../../triangle/mesh.hpp"
#include "../../my_fem/my_fem.hpp"
#include "auxillary.hpp"

/************************************************************************************************************************************/
// Electric field: 	E = -grad(psi)
int compute_ElectricField(	const my_mesh::MeshData &mesh,
				const arma::vec &psi,
				arma::vec &E_x,
				arma::vec &E_y)
{
	arma::vec grad_psi_x, grad_psi_y;
	my_fem::computeGradient(mesh, psi, grad_psi_x, grad_psi_y);
	E_x = -grad_psi_x;
	E_y = -grad_psi_y;
return 0;
}




/************************************************************************************************************************************/
// Diffusion and Convection fluxes

// (1) Diffusion flux: F = - mu * grad(u)
int compute_DiffusionFlux(	const my_mesh::MeshData &mesh,
				const arma::vec &u,
				const arma::vec &mu,			// assumed to be elementwise defined
				arma::vec &F_x,
				arma::vec &F_y)
{
	// 1. Compute gradient of density "u"
	arma::vec grad_u_x, grad_u_y;
	my_fem::computeGradient(mesh, u, grad_u_x, grad_u_y);

	// 2. Assemble coefficient matrix
	arma::mat M = my_fem::assembleMatrixC(mesh, arma::ones<arma::vec>(mesh.num_nodes));

	// 3. Assemble right-hand-side vector
	arma::vec rhs_x(mesh.num_nodes), rhs_y(mesh.num_nodes);
	rhs_x.zeros();	rhs_y.zeros();

	for (int t=0; t<mesh.num_elements; t++)
	{	std::vector<int> neigh_nodes = mesh.topology2to0[t];
		double area_t = mesh.ele_areas[t];
		double mu_t = mu(t);

		for (int i=0; i<neigh_nodes.size(); i++)
		{	int node_index = neigh_nodes[i];
			rhs_x(node_index) += -mu_t*area_t * grad_u_x(node_index) / 3.0;
			rhs_y(node_index) += -mu_t*area_t * grad_u_y(node_index) / 3.0;
		}
	}

	// 4. Solve for flux F = -mu*grad(u)
	F_x = arma::solve(M, rhs_x);
	F_y = arma::solve(M, rhs_y);

return 0;
}


// (2) Convection/Drift flux: F = mu * u * convection
int compute_ConvectionFlux(	const my_mesh::MeshData &mesh,
				const arma::vec &u,				// "u" is the density function
				const arma::vec &conv_x,			// "convection = (conv_x, conv_y)
				const arma::vec &conv_y,
				const arma::vec &mu,				// "mu" is element-wise mobility
				arma::vec &F_x,
				arma::vec &F_y)
{
	// 1. Assemble coefficient matrix corresponding to "int(phi_i, phi_j)"
	arma::mat M = my_fem::assembleMatrixC(mesh, arma::ones<arma::vec>(mesh.num_nodes));

	// 2. Assemble right-hand-side vector
	arma::vec rhs_x(mesh.num_nodes), rhs_y(mesh.num_nodes);
	rhs_x.zeros();	rhs_y.zeros();

	for (int t=0; t<mesh.num_elements; t++)
	{	double area_t = mesh.ele_areas[t];
		double mu_t = mu(t);
		std::vector<int> neigh_nodes = mesh.topology2to0[t];

		for (int i=0; i<neigh_nodes.size(); i++)
		{	int node_index = neigh_nodes[i];
			rhs_x(node_index) += mu_t * area_t * u(node_index) * conv_x(node_index)/3.0;
			rhs_y(node_index) += mu_t * area_t * u(node_index) * conv_y(node_index)/3.0;
		}
	}

	// 3. Solve for solution
	F_x = arma::solve(M, rhs_x);
	F_y = arma::solve(M, rhs_y);
return 0;
}






/************************************************************************************************************************************/
// Electron flux:	Fn = -mu_n(grad(psi)) * grad(n) + mu_n(grad(psi)) * n * grad(psi)
int compute_Flux_n(	const my_mesh::MeshData &mesh,
			const arma::vec &n,
			const arma::vec &psi,
			arma::vec &Fn_x,
			arma::vec &Fn_y)
{
	/** 1. calculate relavant quantities **/
	// 1.1 electric field intensity
	arma::vec E(mesh.num_nodes);
	compute_ElectricFieldAmplitude(	mesh, psi, E);
	// 1.2 gradient of "psi", i.e. "-E" (the convection for electron)
	arma::vec grad_psi_x(mesh.num_nodes), grad_psi_y(mesh.num_nodes);
	my_fem::computeGradient (mesh, psi, grad_psi_x, grad_psi_y);
	// 1.3 elementwise mobility
	arma::vec mu_n(mesh.num_elements);
	compute_MobilityN_elementwise(mesh, E, mu_n);

	/** 2 **/
	// 2.1 Compute diffusion flux
	arma::vec Fn_diff_x, Fn_diff_y;
	compute_DiffusionFlux (mesh, n, mu_n, 					// input
				Fn_diff_x, Fn_diff_y);				// output

	// 2.2 Compute convection flux
	arma::vec Fn_conv_x, Fn_conv_y;
	compute_ConvectionFlux (mesh, n, grad_psi_x, grad_psi_y, mu_n, 		// input
				Fn_conv_x, Fn_conv_y);				// output

	/** 3. Put all together **/
	Fn_x = Fn_diff_x + Fn_conv_x;
	Fn_y = Fn_diff_y + Fn_conv_y;	
 
return 0;
}












/************************************************************************************************************************************/
// Hole flux:		Fp = -mu_p(grad(psi)) * grad(p) + mu_p(grad(psi)) * p * E (E = -grad(psi))
int compute_Flux_p(	const my_mesh::MeshData &mesh,
			const arma::vec &p,
			const arma::vec &psi,
			arma::vec &Fp_x,
			arma::vec &Fp_y)
{
	/** 1. calculate relavant quantities **/
	// 1.1 electric field intensity
	arma::vec E(mesh.num_nodes);
	compute_ElectricFieldAmplitude(	mesh, psi, E);
	// 1.2 Electric field "E" (the convection for hole)
	arma::vec grad_psi_x(mesh.num_nodes), grad_psi_y(mesh.num_nodes);
	my_fem::computeGradient (mesh, psi, grad_psi_x, grad_psi_y);
	arma::vec E_x = - grad_psi_x;
	arma::vec E_y = - grad_psi_y;
	// 1.3 elementwise mobility
	arma::vec mu_p(mesh.num_elements);
	compute_MobilityP_elementwise(mesh, E, mu_p);

	/** 2 **/
	// 2.1 Compute diffusion flux
	arma::vec Fp_diff_x, Fp_diff_y;
	compute_DiffusionFlux (mesh, p, mu_p, 					// input
				Fp_diff_x, Fp_diff_y);				// output

	// 2.2 Compute convection flux
	arma::vec Fp_conv_x, Fp_conv_y;
	compute_ConvectionFlux (mesh, p, E_x, E_y, mu_p, 			// input
				Fp_conv_x, Fp_conv_y);				// output

	/** 3. Put all together **/
	Fp_x = Fp_diff_x + Fp_conv_x;
	Fp_y = Fp_diff_y + Fp_conv_y;	
 
return 0;
}










/************************************************************************************************************************************/
// Exciton flux:	Fx = -mu_x * grad(u)
int compute_Flux_x(	const my_mesh::MeshData &mesh,
			const arma::vec &u,
			arma::vec &Fx_x,					// Fx_x means it's x-component of an exciton flux 
			arma::vec &Fx_y)
{
	// 1. compute relavant quantities
	// 1.1 exciton mobility, element-wise
	arma::vec mu_x;
	compute_MobilityX_elementwise(mesh, mu_x);

	// 2. compute diffusion flux of exciton
	compute_DiffusionFlux(mesh, u, mu_x, Fx_x, Fx_y);

return 0;
}

