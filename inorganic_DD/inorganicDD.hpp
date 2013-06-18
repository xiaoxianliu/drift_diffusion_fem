#ifndef INORGANIC_HPP
#define INORGANIC_HPP

#include <armadillo>
#include "../triangle/mesh.hpp"

// Main loop of Gummel's iteration
int gummelIteration(	const my_mesh::MeshData &mesh,
			arma::vec &psi,
			arma::vec &n,
			arma::vec &p);


// Solvers for each equations
arma::vec solveNContinuityEq(	const my_mesh::MeshData &mesh,
				const arma::vec &input_psi,
				const arma::vec &input_n,
				const arma::vec &input_p,
				const arma::vec &n_D);				// Dirichlet boundary condition
arma::vec solvePContinuityEq(	const my_mesh::MeshData &mesh,
				const arma::vec &input_psi,
				const arma::vec &input_n,
				const arma::vec &input_p,
				const arma::vec &p_D);				// Dirichlet boundary condition

arma::vec solveNonlinearPoissonEq(	const my_mesh::MeshData &mesh,
					const arma::vec &input_psi,
					const arma::vec &input_n,
					const arma::vec &input_p,
					const arma::vec &psi_D,
					const arma::vec &C	);			// "C" is arma::vec form of doping density










/** Physical settings ****************************************************/
arma::vec doping_Vec	(const my_mesh::MeshData &mesh, double C_p=-1.0, double C_n=1.0);	// Dopant density
arma::vec builtinPsi_Vec(const arma::vec &C);							// Built-in potential


/** Functions computing reaction rates ***********************************/
arma::vec recombSRH_Vec(const arma::vec n, const arma::vec p);
double generation_Func(double x, double y);


/** Dirichlet boundary conditions ****************************************/
arma::vec n_D_Vec(arma::vec C);
arma::vec p_D_Vec(arma::vec C);
arma::vec psi_D_Vec(	const my_mesh::MeshData &mesh, 
			const arma::vec &psi_bi,
			double psi_applied);

#endif

