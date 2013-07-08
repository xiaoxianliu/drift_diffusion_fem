#ifndef AUXILLARY_HPP
#define AUXILLARY_HPP


#include <cmath>
#include <vector>

#include <armadillo>

#include "../../triangle/mesh.hpp"
#include "../../my_fem/my_fem.hpp"

#include "parameters.hpp"


/*******************************************************************************************************************/
// Compute magnitude of electrical field |E|
int compute_ElectricFieldAmplitude(	const my_mesh::MeshData &mesh,
					const arma::vec &psi,
					arma::vec &E);

/*******************************************************************************************************************/
// mobility of electrons and holes
// 1. element-wise mobility
int compute_MobilityN_elementwise (	const my_mesh::MeshData &mesh,
					const arma::vec &E,				// electric field intensity
					arma::vec &mu_n_elementwise);

int compute_MobilityP_elementwise (	const my_mesh::MeshData &mesh, 
					const arma::vec &E, 				// electric field intensity
					arma::vec &mu_p_elementwise);
int compute_MobilityX_elementwise (	const my_mesh::MeshData &mesh,
					arma::vec &mu_x_elementwise);


/*******************************************************************************************************************/
// Compute "outward" normal derivative of given vector "u" on interface w.r.t. either left region (marked "1") or right region (marked "2")
int compute_dpsi_dnu_interface(	const my_mesh::MeshData &mesh,
					const arma::vec &psi,
					const arma::vec &n,
					const arma::vec &p,
					const arma::vec &u,
					arma::vec &dpsi_dnu_1,
					arma::vec &dpsi_dnu_2);
/*******************************************************************************************************************/
// Compute surface rate of exciton dissociation on interface; (surface rate = volume rate * interface_width)
// it's a vector taking non-zero values only on the interface
int compute_ExcitonDissociationRate(	const my_mesh::MeshData &mesh,
					const arma::vec &psi,
					const arma::vec &n,
					const arma::vec &p,
					const arma::vec &u,
					arma::vec &k_diss_interface);

/*******************************************************************************************************************/
// Compute surface rate of recombination on interface; (surface rate = volume rate * interface_width)
// it's vector taking non-zero values only on the interface
int compute_RecombinationRate(	const my_mesh::MeshData &mesh,
				const arma::vec &E,				// electric field amplitude
				arma::vec &recomb_interface);





/*************************************************************************************************************************/
// Compute arma::vec of "photo generation function"
int compute_PhotoGenerationVec(	const my_mesh::MeshData &mesh, arma::vec &Q_vec);






#endif
