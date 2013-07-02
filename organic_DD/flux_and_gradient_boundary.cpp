#include <iostream>
#include <cstdlib>

#include <armadillo>

#include "../triangle/mesh.hpp"
#include "../my_fem/my_fem.hpp"

/************************************************************************************************************************************/
// Electric field: 	dE/dnu = -dpsi/dnu
// (1) On Gamma_anode:	dE/dnu1
// (2) On Gamma_interface:	dE/dnu1, dE/dnu2
// (3) On Gamma_cathod:	dE/dnu2

/************************************************************************************************************************************/
// Electron flux:	Fn_dot_nu
// (1) On Gamma_anode:	Fn_dot_nu1
int compute_BoundaryFlux_n(	const my_mesh::MeshData &mesh,
				const arma::vec &n,
				const arma::vec &psi,
				int contact_marker,
				arma::vec &Fn_nu)
{
	// 0. Check if "contact_marker" is either 1 or 3; if not, exit with error message
	if (contact_marker != 1 && contact_marker != 3)
	{	std::cout << "\"Contact_marker\" has to be one of \"1\" and \"3\"; \"1\" for anode and \"3\" for cathode. Exiting.\n";
		exit(1);
	}
	// 1. Assemble coefficient matrix on boundary marked by "contact_marker"

return 0;
}
// (2) On Gamma_interface:	Fn_dot_nu1, Fn_dot_nu2
// (3) On Gamma_cathode:	Fn_dot_nu2


/************************************************************************************************************************************/
// Hole flux:		Fp_dot_nu
// (1) On Gamma_anode:	Fp_dot_nu1
// (2) On Gamma_interface:	Fp_dot_nu1, Fp_dot_nu2
// (3) On Gamma_cathode:	Fp_dot_nu2

/************************************************************************************************************************************/
// Exciton flux:	Fx_dot_nu
// (1) On Gamma_anode:	Fx_dot_nu1
// (2) On Gamma_interface:	Fx_dot_nu1, Fx_dot_nu2
// (3) On Gamma_cathode:	Fx_dot_nu2
