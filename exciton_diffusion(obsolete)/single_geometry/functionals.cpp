#include <armadillo>
#include "mesh.hpp"
#include "exciton_diffusion.hpp"

double interfaceIntegral(const MeshData& mesh, const arma::vec& u, const arma::vec& d)
{	double sum=0;;

	for (int e=0; e<mesh.num_edges; e++)
	{
		if (mesh.edge_markers[e]==5)
		{	int v0 = mesh.edges[2*e];
			int v1 = mesh.edges[2*e+1];
			double length = mesh.edge_lengths[e];

			sum += length * (u(v0)*d(v0) + u(v1)*d(v1)) / 2.0;		// Trapezoidal rule for 1st order approx.
		}
	}
	return sum;
}
