#ifndef MAIN_HPP
#define MAIN_HPP

#include "../triangle/mesh.hpp"
#include <armadillo>

int solve_GummelIteration(	const my_mesh::MeshData &mesh,
				arma::vec &psi,
				arma::vec &n,
				arma::vec &p,
				arma::vec &u,
				double applied_psi = 0.0			// applied potential at "anode"
			);


#endif
