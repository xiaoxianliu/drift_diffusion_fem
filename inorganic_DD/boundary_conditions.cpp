#include <cmath>
#include <armadillo>
#include "../triangle/mesh.hpp"
#include "parameters.hpp"

// Dirichlet boundary conditions for n, p, and psi
arma::vec n_D_Vec(arma::vec C)
{using parameters::delta_squared;

	arma::vec n_D(C.n_cols);
	for (int i=0; i<C.n_cols; i++)
	{	double c = C(i);
		n_D(i) = 0.5*(c + sqrt(c*c + 4*delta_squared*delta_squared));
	}
	return n_D;
}






arma::vec p_D_Vec(arma::vec C)
{using parameters::delta_squared;

	arma::vec p_D(C.n_cols);
	for (int i=0; i<C.n_cols; i++)
	{	double c = C(i);
		p_D(i) = 0.5*(-c+ sqrt(c*c + 4*delta_squared*delta_squared));
	}
	return p_D;
}







arma::vec psi_D_Vec(	const my_mesh::MeshData &mesh, 
			const arma::vec &psi_bi,
			double psi_applied)
{
using parameters::delta_squared;

	arma::vec psi_D;
	psi_D = psi_bi;
	for (int i=0; i<mesh.num_nodes; i++)
	{	std::vector<int> neigh_eles = mesh.topology0to2[i];
		for (int j=0; j<neigh_eles.size(); j++)
		{	int t = neigh_eles[j];
			if (mesh.element_markers[t]==1)
			{	psi_D(i) += psi_applied;
				break;
			}
		}
	}

	return psi_D;
}






