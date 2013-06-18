#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <armadillo>

#include "../triangle/mesh.hpp"
#include "parameters.hpp"

// 1. Dopant density
arma::vec doping_Vec(	const my_mesh::MeshData &mesh, double C_p, double C_n)
{
	int N = mesh.num_nodes;
	arma::vec C(N);

	for (int i=0; i<N; i++)					// loop through all vertices
	{
		// (1) Find the correct marker
		std::vector<int> neigh_ele = mesh.topology0to2[i];
		int marker=1;
		for (int j=0; j<neigh_ele.size(); j++)
		{	int t = neigh_ele[j];
			if (mesh.element_markers[t]!=1)
			{	marker = mesh.element_markers[t];
				if (marker!=2)
				{	std::cout << "Marker of the "<<t<<"-th element is "<<marker
						<<"; it should be either 1 or 2\n";
					exit(1);
				}
				break;
			}
		}
		// (2) Set the doping density
		if (marker == 1)	C(i) = C_p;
		else	C(i) = C_n;
	}
	return C;
}



// 2. Built-in potential
arma::vec builtinPsi_Vec(const arma::vec &C)		// given doping density "C" in the form of arma::vec
{
using parameters::delta_squared;

	arma::vec psi_bi(C.n_rows);
	for (int i=0; i<C.n_rows; i++)
	{	psi_bi(i) = log( (C(i) + sqrt(C(i)*C(i) + 4*delta_squared*delta_squared) )
				/(2*delta_squared) );
	}
	return psi_bi;
}	

