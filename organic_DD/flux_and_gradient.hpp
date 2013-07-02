#ifndef FLUX_GRAD_HPP
#define FLUX_GRAD_HPP


#include <armadillo>

#include "../triangle/mesh.hpp"
//#include "../my_fem/my_fem.hpp"
//#include "auxillary.hpp"

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





#endif
