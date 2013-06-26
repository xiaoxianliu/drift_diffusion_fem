#include <cmath>
#include <vector>

#include <armadillo>

#include "../triangle/mesh.hpp"
#include "../my_fem/my_fem.hpp"
#include "parameters.hpp"
//#include "main.hpp"

/*******************************************************************************************************************/
// Compute magnitude of electrical field |E|
int compute_ElectricFieldAmplitude(	const my_mesh::MeshData &mesh,
					const arma::vec &psi,
					arma::vec &E);

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
