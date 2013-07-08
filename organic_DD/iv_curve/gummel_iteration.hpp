#ifndef GUMMEL_HPP
#define GUMMEL_HPP

#include <armadillo>
#include "../../triangle/mesh.hpp"


// 1. Declaration of individual solvers for each equation
// 1.1 continuity equation for electron
int solve_ContinuityEq_n(const my_mesh::MeshData &mesh, 
			const arma::vec &input_psi,
			const arma::vec &input_n,
			const arma::vec &input_p,
			const arma::vec &input_u,
			arma::vec &output_n);

// 1.2 continuity equation for holes
int solve_ContinuityEq_p(const my_mesh::MeshData &mesh, 
			const arma::vec &input_psi,
			const arma::vec &input_n,
			const arma::vec &input_p,
			const arma::vec &input_u,
			arma::vec &output_p);

// 1.3 continuity equation for excitons
int solve_ContinuityEq_x(	const my_mesh::MeshData &mesh,
				const arma::vec &input_psi,
				const arma::vec &input_n,
				const arma::vec &input_p,
				const arma::vec &input_u,
				arma::vec &output_u);
// 1.4 nonlinear poisson's equation
int solve_NonlinearPoissonEq(	const my_mesh::MeshData &mesh,
				const arma::vec &input_psi,
				const arma::vec &input_n,
				const arma::vec &input_p,
				const arma::vec &input_u,
				arma::vec &output_psi,
				double applied_psi);


// 2. Declare the function implementing Gummel's iteration
int solve_GummelIteration(	const my_mesh::MeshData &mesh,
				arma::vec &psi,
				arma::vec &n,
				arma::vec &p,
				arma::vec &u,
				double applied_psi			// applied potential at "anode"
			);


#endif
