#include <cmath>

double f(double x, double y)
{	double sigma = 0.1;
	return 1 + exp( -0.5*pow(x,2)/pow(sigma,2) - 0.5*pow(y-0.5,2)/pow(sigma,2) );
}

double grad_f_x(double x, double y)
{	double result, sigma=0.1;
	result = exp( -0.5*pow(x,2)/pow(sigma,2) - 0.5*pow(y-0.5,2)/pow(sigma,2) )
		*(-1/pow(sigma,2)) * x;
	return result;
}

double grad_f_y(double x, double y)
{	double result, sigma=0.1;
	result = exp( -0.5*pow(x,2)/pow(sigma,2) - 0.5*pow(y-0.5,2)/pow(sigma,2) )
		*(-1/pow(sigma,2)) * (y - 0.5);
	return result;
}
