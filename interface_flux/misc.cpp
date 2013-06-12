/*** Functions used only locally ****/





/****** Functions related to state equation of unknown "u" ******/

/* FEM formulation of state equation:
/* 	int_Omega ( a * <grad(u), grad(v)> + c * u*v) dx + int_Gamma (d * u*v ) ds = int_Omega ( g * v) dx  
/* where
/*	Omega: entire domain = union(Omega1, Omega2)
/*	Gamma: interface between Omega1 and Omega2
/*	v: test function
*/

double func_a(double x, double y)				// defusion coefficient
{	return 1.0;
}
double func_c(double x, double y)				// body decay rate
{	return 1.0;
}
double func_d(double x, double y)				// interface decay rate
{	return 1.0;
}
double func_g(double x, double y)				// body source
{	return 1.0;
}
double uD(double x, double y)					// Dirichlet boundary condition
{	return 0.0;
}








