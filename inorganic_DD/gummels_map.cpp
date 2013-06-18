#include <iostream>
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
	double psi_applied = 1.0;						// applied potential on p-end
	arma::vec psi_bi = builtinPsi_Vec(C);					// built-in potential


	/* 2. Compute Dirichlet boundary conditions */
	arma::vec n_D = n_D_Vec(C);
	arma::vec p_D = p_D_Vec(C);
	arma::vec psi_D = psi_D_Vec(mesh, psi_bi, psi_applied);


/*	psi = psi_D;	n = n_D;	p = p_D;
	plot_ArmaVec(mesh, C, "dopant_density");
	plot_ArmaVec(mesh, psi_bi, "built-in_potential");
	plot_ArmaVec(mesh, psi, "initial_psi");
	plot_ArmaVec(mesh, n, "initial_n");
	plot_ArmaVec(mesh, p, "initial_p");
*/



	/* 3. Gummel's iteration */
	arma::vec prev_psi = psi_D, prev_n = n_D, prev_p = p_D;			// initialization
	
	arma::vec 	new_n(mesh.num_nodes), 
			new_p(mesh.num_nodes), 
			new_psi(mesh.num_nodes);

	int max_iter = 5;
	double gummel_err = 1.0, gummel_tol = 1e-3;
	for (int i=0; i<max_iter; i++)
	{

		std::cout << "Gummel iteration " << i << ":\n";




		// 3.1 compute solution of N-continuity equation 
		new_n = solveNContinuityEq (mesh, prev_psi, prev_n, prev_p, n_D);

		// 3.2 compute solution of P-continuity equation 
		new_p = solvePContinuityEq (mesh, prev_psi, prev_n, prev_p, p_D);

		// 3.3 compute solution of (nonlinear) poisson's equation 
		new_psi = solveNonlinearPoissonEq (mesh, prev_psi, new_n, new_p, psi_D, C);

		// 3.4 Compute Gummel's error 
		double max_dn = arma::norm( arma::abs(new_n - prev_n), "inf");
		double max_dp = arma::norm( arma::abs(new_p - prev_p), "inf");
		double max_dpsi = arma::norm( arma::abs(new_psi - prev_psi), "inf");
		gummel_err = fmax( fmax(max_dn, max_dp), max_dpsi);

		// 3.5 Replace old vectors by the newly computed
		prev_n = new_n;	prev_p = new_p;	prev_psi = new_psi;


		std::cout << "\t gummel_err = " << gummel_err << "\n";

		if (gummel_err < gummel_tol)
			break;
	}


	psi = new_psi;	n = new_n;	p = new_p;

	plot_ArmaVec(mesh, psi, "psi");
	plot_ArmaVec(mesh, n, "n");
	plot_ArmaVec(mesh, p, "p");

	
	return 0;
}
