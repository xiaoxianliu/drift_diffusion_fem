#include <armadillo>
#include "fem_assemble.hpp"
#include "../triangle/mesh.hpp"

namespace linear_fem
{
arma::vec interpolateFunction(const my_mesh::MeshData &mesh, double (*f)(double, double))
{	int N = mesh.num_nodes;

	arma::vec u(N);

	for (int i=0; i<N; i++)
	{	double x = mesh.nodes[i][0];
		double y = mesh.nodes[i][1];
		u(i) = f(x,y);
	}
	return u;
}

}
