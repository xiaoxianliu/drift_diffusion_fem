#include <cmath>

namespace parameters
{

// References:Barker, Ramsdale, Greenham, 2003 
// Buxton and Clarke, 2007
// Falco, Porro, Sacco, and Verri, 2012


/**** Universal constant (not supposed to be changed)**/
double epsilon_vacuum = 8.854187817e-12;		// vacuum permittivity, in "F_m^(-1)"
double k_boltzmann = 1.3806488e-23;			// Boltzmann constant, in "m2_kg_s^(-2)_K^(-1)
double q_unit = 1.602176565e-19;			// unit charge, in "C"




// Temperature and thermal potential 
double temperature = 300;				// Kelvin
double U_T = k_boltzmann * temperature / q_unit;	// thermal potential, in "V"; should be around 0.026V



/******** typical values of different physical quantities ********/
// To de-dimensionlize quantities below
double typical_length = 2e-7;				// device dimension is assumed to be 200 nm = 2e-7m
double typical_mobility = 1e-9;				// measure of mobility is 1e-9 m^2_V^(-1)_s^(-1)
double typical_density = 1e20;				// density of particles: 1e20 m^(-3)
double typical_time = pow( typical_length, 2) / (typical_mobility * U_T);
 							// time scale associated with diffusion/drift process



/******* device geometry *****************************************************************************/
// Physical value
double interface_width = 2e-9;				// 2 nm width for the interface region
// Dimensionless value
double h_interface = interface_width / typical_length;






/**** (relative) permittivity ************************************************************************/
double epsilon_acceptor = 4.0;				// relative permittivity "epsilon_r"
double epsilon_donor = 4.0;
double epsilon_rel = 0.5*(epsilon_acceptor + epsilon_donor);


/**** Exciton characteristics ************************************************************************/
// Dimensional values
double typical_exciton_radius = 1e-9;			// typical exciton radius is 1nm
double exciton_lifetime = 1e-6;				// s
// Dimensionless values
double a_x = typical_exciton_radius / typical_length;	// exciton size
double tau_x = exciton_lifetime / typical_time;		// exciton lifetime



/**** Mobility (diffusivity) etc ********************************************************************/
// 1. Physical values
// 1.1 Exciton
double exciton_diffusivity = 1e-10;			// m^2_s^(-1);	
double exciton_mobility = exciton_diffusivity / U_T;
double exciton_diffusion_length = sqrt(exciton_diffusivity * exciton_lifetime);		// about 10nm
							// around 10nm = 1e-8m
// 1.2 Electron
double electron_mobility_donor = 1e-9;			// Zero-field mobilities, in "m^2_V^(-1)_s^(-1)"
double electron_mobility_acceptor = 1e-8;	
double electron_mobility_field_dependent_coefficient = 1.55e-3;		// m^(1/2)_V^(-1/2); "mu_n = mu0_n exp(gamma_n * sqrt(|E|))"

// 1.3 Hole
double hole_mobility_donor = 2e-8;
double hole_mobility_acceptor = 2e-9;
double hole_mobility_field_dependent_coefficient = 3e-4;		// m^(1/2)_V^(-1/2); "mu_p = mu0_p exp(gamma_p * sqrt(|E|))"


// 2. Dimensionless values
double mu_x = exciton_mobility / typical_mobility;		// exciton mobility
double mu_n_donor = electron_mobility_donor / typical_mobility;	// electron mobility in donor and acceptor
double mu_n_acceptor = electron_mobility_acceptor / typical_mobility;
double mu_p_donor = hole_mobility_donor / typical_mobility;	// hole mobility in donor and acceptor
double mu_p_acceptor = hole_mobility_acceptor / typical_mobility;
double gamma_n = electron_mobility_field_dependent_coefficient * sqrt( U_T / typical_length );	
								// field dependent parameter for electron mobilities
double gamma_p = hole_mobility_field_dependent_coefficient * sqrt( U_T / typical_length );	
								// field dependent parameter for hole mobilities






/**** Reactions *********************************************************************************/
// 1. Physical values
// 1.1 Light absorption
double photon_absorption_rate = 1e25;			// m^(-3)_s^(-1), i.e. counts per unit volume per second
							// Ref: Falco, Porro, etc
double absorption_coefficient = 3e6;			// "m^(-1)"; Q ~ exp(-alpha*x)
double absorption_length = 1.0/absorption_coefficient;	// "m"
// 1.2 Recombination of electron-hole pair
double typical_recombination_rate = q_unit * typical_mobility / epsilon_vacuum ;
// 1.3 Exciton dissociation
double exciton_dissociation_rate_zero_field = 1e8;			// s^(-1)

// 2. Dimensionless form
// Light absorption
double Q0 = photon_absorption_rate / (typical_density / typical_time);
double alpha = absorption_coefficient * typical_length;		// In dimensionless form, Q ~ Q0*exp(- alpha*x)
// recombination
double recomb_coefficient = typical_recombination_rate * (typical_density * typical_time);
// exciton dissociation
double k_diss_0 = exciton_dissociation_rate_zero_field * typical_time;		// zero-field dissociation







/**** Dimensionless quantities ****/
double lambda_squared = epsilon_vacuum * U_T / (q_unit * typical_density * pow(typical_length, 2));




/***** Dirichlet boundary values ************************************************************************/
// 1. Physical values
// 1.1 electrical potential
double builtin_potential_anode = -0.4;
double builtin_potential_cathode = 0.0;
double builtin_potential = builtin_potential_anode - builtin_potential_cathode;

// 1.2 electron
double electron_density_anode = 0.0;
double electron_density_cathode = 1e22;
// 1.3 hole
double hole_density_anode = 1e22;
double hole_density_cathode = 0.0;
// 1.4 exciton
double exciton_density_anode = 0.0;
double exciton_density_cathode = 0.0;

// 2. Dimensionless values
double psi_A = builtin_potential_anode / U_T;
double psi_C = builtin_potential_cathode / U_T;
double psi_bi = psi_A - psi_C;

double n_A = electron_density_anode / typical_density;
double n_C = electron_density_cathode / typical_density;

double p_A = hole_density_anode / typical_density;
double p_C = hole_density_cathode / typical_density;

double u_A = exciton_density_anode / typical_density;
double u_C = exciton_density_cathode / typical_density;





}
