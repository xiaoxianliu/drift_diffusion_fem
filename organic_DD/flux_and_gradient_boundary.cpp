#include <iostream>

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
