#include <cmath>
#include <armadillo>
#include "../triangle/mesh.hpp"
#include "inorganicDD.hpp"

int gummelIteration(	const my_mesh::MeshData &mesh,
			arma::vec &psi,
			arma::vec &n,
			arma::vec &p)
{
	/* 1. Initialize physical input and unknowns */
	double C_p = -1.0, C_n = 1.0;
	arma::vec C = doping_Vec(mesh, C_p, C_n);				// Dopant density
	double psi_applied = 0.0;						// applied potential on p-end
	arma::vec psi_bi = builtinPsi_Vec(C);					// built-in potential

	/* 2. Compute Dirichlet boundary conditions */
	arma::vec n_D = n_D_Vec(C);
	arma::vec p_D = p_D_Vec(C);
	arma::vec psi_D = psi_D_Vec(mesh, psi_bi, psi_applied);
/*
	psi = psi_D;	n = n_D;	p = p_D;
*/
	/* 3. Gummel's iteration */
/*	arma::vec prev_psi = psi_D, prev_n = n_D, prev_p = p_D;			// initialization
	
	arma::vec new_n, new_p, new_psi;
	int max_iter = 10;
	double gummel_err = 1.0, gummel_tol = 1e-3;
	for (int i=0; i<max_iter; i++)
	{
		// 3.1 compute solution of N-continuity equation 
		arma::vec new_n = solveNContinuityEq (mesh, prev_psi, prev_n, prev_p, n_D);

		// 3.2 compute solution of P-continuity equation 
		arma::vec new_p = solvePContinuityEq (mesh, prev_psi, prev_n, prev_p, p_D);

		// 3.3 compute solution of (nonlinear) poisson's equation 
		arma::vec new_psi = solveNonlinearPoissonEq (mesh, prev_psi, new_n, new_p, psi_D, C);

		// 3.4 Compute Gummel's error 
		double max_dn = arma::norm( arma::abs(new_n - prev_n), "inf");
		double max_dp = arma::norm( arma::abs(new_p - prev_p), "inf");
		double max_dpsi = arma::norm( arma::abs(new_psi - prev_psi), "inf");
		gummel_err = fmax( fmax(max_dn, max_dp), max_dpsi);

		if (gummel_err < gummel_tol)
			break;
	}


	psi = new_psi;	n = new_n;	p = new_p;
*/
	return 0;
}
