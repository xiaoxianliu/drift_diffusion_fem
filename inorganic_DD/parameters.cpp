#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <cmath>

namespace parameters
{
// Prescribed physical parameters
double n_intrinsic = 1E10;							// 3-D intrinsic carrier density in cm^(-3)
double N_typical = 1E15;							// Typical doping level in cm^(-3)
double l_typical = 2E-3;							// Typical length of devices in cm. 
double mu_typical = 1E3;    							// Typical mobility, in cm^2/(Volt*sec)
double epsilon_0 = 8.854E-14;							// Vacuum permittivity in F/cm
double epsilon_relative = 11.68;						// Relative permittivity of Silicon
double q_unit = 1.602E-19;							// Unit charge of electron (absolute value) in Coulomb
double temperature = 300;							// Room temperature in Kelvin
double k_boltzman = 1.38E-23;							// Boltzmann constant in ( m^2 * kg * s^(-2) * K^(-1))



// Computed parameters
double epsilon = epsilon_0 * epsilon_relative;	// Permittivity of Silicon in physical dimension, same dimension as "epsilon_0"
double U_T = k_boltzman * temperature / q_unit;		// Thermal potential in "kT/q" in Volt
double l_debye = sqrt(epsilon*U_T/(q_unit * N_typical));	// Debye length, in cm
double tau_typical = (l_typical*l_typical) / (mu_typical * U_T);		// Time unit, in "sec"
double D_typical = U_T * mu_typical;				// The typical value for diffusion coefficient.


// Computed dimensionless parameters
double lambda_squared = pow(l_debye / l_typical, 2);
double delta_squared = n_intrinsic / N_typical;


// The quantities below are only present in the continuity equations, not the nonlinear Poisson equations
// dimensionless mobilities
double mu_n = 1.345;					// in "mu_tilt" defined above.
double mu_p = 0.477;


// dimensionless minority carriers' lifetime
double tau_n = 1000.0;					// in "tau_tilt" computed above
double tau_p = 1000.0;

// dimensionless generation rate
//double generation = 0.0;





}



#endif
