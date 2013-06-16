
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






/******* Functions related to adjoint equation of unknown "xi" *********/
/* FEM formaulation of adjoint equation:
/* 	int_Omega ( a_adjoint * <grad(u), grad(v)> + c_adjoint * u*v) dx + int_Gamma (d_adjoint * u*v ) ds = int_Gamma ( 1 * v) dx  
/* where
/*	Omega: entire domain = union(Omega1, Omega2)
/*	Gamma: interface between Omega1 and Omega2
/*	v: test function
/* Note: the interface source term on right-hand side
*/


double func_a_adjoint(double x, double y)
{	return 1.0;
}
double func_c_adjoint(double x, double y)
{	return 1.0;
}
double func_d_adjoint(double x, double y)
{	return 1.0;
}
double xiD(double x, double y)				// xi_D is Dirichlet boundary condition
{	return 0.0;
}


