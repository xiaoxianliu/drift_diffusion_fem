#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <cmath>

namespace parameters
{
// Prescribed physical parameters
extern double n_intrinsic;							// 3-D intrinsic carrier density in cm^(-3)
extern double N_typical;							// Typical doping level in cm^(-3)
extern double l_typical;							// Typical length of devices in cm. 
extern double mu_typical;    							// Typical mobility, in cm^2/(Volt*sec)
extern double epsilon_0;							// Vacuum permittivity in F/cm
extern double epsilon_relative;						// Relative permittivity of Silicon
extern double q_unit;							// Unit charge of electron (absolute value) in Coulomb
extern double temperature;							// Room temperature in Kelvin
extern double k_boltzman;							// Boltzmann constant in ( m^2 * kg * s^(-2) * K^(-1))



// Computed parameters
extern double epsilon;			// Permittivity of Silicon in physical dimension, same dimension as "epsilon_0"
extern double U_T;			// Thermal potential in "kT/q" in Volt
extern double l_debye;			// Debye length, in cm
extern double tau_typical;		// Time unit, in "sec"
extern double D_typical;		// The typical value for diffusion coefficient.


// Computed dimensionless parameters
extern double lambda_squared;
extern double delta_squared;


// The quantities below are only present in the continuity equations, not the nonlinear Poisson equations
// dimensionless mobilities
extern double mu_n;					// in "mu_tilt" defined above.
extern double mu_p;


// dimensionless minority carriers' lifetime
extern double tau_n;					// in "tau_tilt" computed above
extern double tau_p;





}



#endif
