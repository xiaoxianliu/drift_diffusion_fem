#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <cmath>
namespace parameters
{

// References:Barker, Ramsdale, Greenham, 2003 
// Buxton and Clarke, 2007
// Falco, Porro, Sacco, and Verri, 2012


/**** Universal constant (not supposed to be changed)**/
extern double epsilon_vacuum;					// vacuum permittivity, in "F_m^(-1)"
extern double k_boltzmann;					// Boltzmann constant, in "m2_kg_s^(-2)_K^(-1)
extern double q_unit;						// unit charge, in "C"



// Temperature and thermal potential 
extern double temperature;						// Kelvin
extern double U_T;							// thermal potential, in "V"; should be around 0.026V




/******* device geometry *****************************************************************************/
// Dimensionless value of interface width
extern double h_interface;




/**** (relative) permittivity ***********************/
extern double epsilon_acceptor;						// relative permittivity "epsilon_r"
extern double epsilon_donor;
extern double epsilon_rel;						// average relative permittivity "epsilon_r"


/**** Exciton characteristics *****/

// Dimensionless values
extern double a_x;	// exciton size
extern double tau_x;		// exciton lifetime



/**** Mobility (diffusivity) etc ***********/
// 2. Dimensionless values
extern double mu_x;		// exciton mobility
extern double mu_x_donor;
extern double mu_x_acceptor;
extern double mu_n_donor;	// electron mobility in donor and acceptor
extern double mu_n_acceptor;
extern double mu_p_donor;	// hole mobility in donor and acceptor
extern double mu_p_acceptor;
extern double gamma_n;						// field dependent parameter for electron mobilities
extern double gamma_p;						// field dependent parameter for hole mobilities






/**** Reactions ****************/
// 2. Dimensionless form
// Light absorption
extern double Q0;
extern double alpha;		// In dimensionless form, Q ~ Q0*exp(- alpha*x)
// recombination
extern double recomb_coefficient;				// recomb_rate = recomb_coefficient * (mu_n + mu_p) / <epsilon_rel>
// exciton dissociation
extern double k_diss_0;







/**** Dimensionless quantities ****/
extern double lambda_squared;






/***** Dirichlet boundary values *********/
// 2. Dimensionless values
extern double psi_A;
extern double psi_C;
extern double psi_bi;

extern double n_A;
extern double n_C;

extern double p_A;
extern double p_C;

extern double u_A;
extern double u_C;



}

#endif
