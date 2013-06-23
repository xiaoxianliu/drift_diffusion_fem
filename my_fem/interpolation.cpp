#include <armadillo>
#include "my_fem.hpp"
#include "../triangle/mesh.hpp"

namespace my_fem
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


arma::vec interpolateConstant(const my_mesh::MeshData &mesh, double C)
{
	int N = mesh.num_nodes;
	arma::vec u(N);
	for (int i=0; i<N; i++)
	{	u(i) = C;
	}
	return u;
}

}
