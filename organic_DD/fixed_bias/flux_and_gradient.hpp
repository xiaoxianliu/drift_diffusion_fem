#ifndef FLUX_GRAD_HPP
#define FLUX_GRAD_HPP


#include <armadillo>

#include "../../triangle/mesh.hpp"


/************************************************************************************************************************************/
// Electric field: 	E = -grad(psi)
int compute_ElectricField(	const my_mesh::MeshData &mesh,
				const arma::vec &psi,
				arma::vec &E_x,
				arma::vec &E_y);



/************************************************************************************************************************************/
// Diffusion and Convection fluxes

// (1) Diffusion flux: F = - mu * grad(u)
int compute_DiffusionFlux(	const my_mesh::MeshData &mesh,
				const arma::vec &u,
				const arma::vec &mu,			// assumed to be elementwise defined
				arma::vec &F_x,
				arma::vec &F_y);


// (2) Convection/Drift flux: F = mu * u * convection
int compute_ConvectionFlux(	const my_mesh::MeshData &mesh,
				const arma::vec &u,				// "u" is the density function
				const arma::vec &conv_x,			// "convection = (conv_x, conv_y)
				const arma::vec &conv_y,
				const arma::vec &mu,				// "mu" is element-wise mobility
				arma::vec &F_x,
				arma::vec &F_y);

/************************************************************************************************************************************/
// Electron flux:	Fn = -mu_n(grad(psi)) * grad(n) + mu_n(grad(psi)) * n * grad(psi)
int compute_Flux_n(	const my_mesh::MeshData &mesh,
			const arma::vec &n,
			const arma::vec &psi,
			arma::vec &Fn_x,
			arma::vec &Fn_y);

/************************************************************************************************************************************/
// Hole flux:		Fp = -mu_p(grad(psi)) ( grad(p) + p*grad(psi))
int compute_Flux_p(	const my_mesh::MeshData &mesh,
			const arma::vec &p,
			const arma::vec &psi,
			arma::vec &Fp_x,
			arma::vec &Fp_y);

/************************************************************************************************************************************/
// Exciton flux:	Fx = -mu_x * grad(u)
int compute_Flux_x(	const my_mesh::MeshData &mesh,
			const arma::vec &u,
			arma::vec &Fx_x,					// Fx_x means it's x-component of an exciton flux 
			arma::vec &Fx_y);




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
				arma::vec &Fn_nu1);

// (2) On Gamma_interface:	Fn_dot_nu1, Fn_dot_nu2
// (3) On Gamma_cathode:	Fn_dot_nu2

/************************************************************************************************************************************/
// Normal component of hole flux:		Fp_dot_nu
// (1) On Gamma_anode:	Fp_dot_nu1
int compute_Boundary1Flux_p(	const my_mesh::MeshData &mesh,
				const arma::vec &p,
				const arma::vec &psi,
				arma::vec &Fp_nu1);
// (2) On Gamma_interface:	Fp_dot_nu1, Fp_dot_nu2
// (3) On Gamma_cathode:	Fp_dot_nu2

/************************************************************************************************************************************/
// Normal component of exciton flux:		Fx_dot_nu
// (1) On Gamma_anode:	Fx_dot_nu1
// (2) On Gamma_interface:	Fx_dot_nu1, Fx_dot_nu2
// (3) On Gamma_cathode:	Fx_dot_nu2










/************************************************************************************************************************************/
// averaged current density
// (1) electron current on anode
double compute_Boundary1CurrentDensity_n(	const my_mesh::MeshData &mesh,
						const arma::vec &n,
						const arma::vec &psi);
// (2) hole current on anode
double compute_Boundary1CurrentDensity_p(	const my_mesh::MeshData &mesh,
						const arma::vec &p,
						const arma::vec &psi);
#endif
