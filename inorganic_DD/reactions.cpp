#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <armadillo>
#include "../triangle/mesh.hpp"
#include "parameters.hpp"

// 1. SRH recombination rate
arma::vec recombSRH_Vec(	const arma::vec n,
				const arma::vec p)
{
using parameters::delta_squared;
using parameters::tau_n;
using parameters::tau_p;

	int N = n.n_rows;
	if ( N!=p.n_rows )
	{	std::cout << "In the function computing SRH recombination rate...\n";
		std::cout << "Dimension of \"n\" is " << n.n_rows << " dimension of \"p\"  is " << p.n_rows << "\n";
		std::cout << "They have to be equal! Exiting... \n";
		exit(1);
	}
	// Compute element-wise value for recombination
	arma::vec R(N);	R.zeros();
	for (int i=0; i<N; i++)
	{	double numerator = n(i) * p(i) - pow(delta_squared, 2);
		double denominator = tau_p * (n(i) + delta_squared) + tau_n * (p(i) + delta_squared);
		R(i) = numerator/denominator;
	}
	return R;
}


// 2. generation rate
double generation_Func(double x, double y)
{
	return 0.0;
}



