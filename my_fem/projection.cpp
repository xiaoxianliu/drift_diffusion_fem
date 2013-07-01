#include <vector>
#include <armadillo>
#include "../triangle/mesh.hpp"

namespace my_fem
{
/**** Compute int_Omega (u * phi_i)dx for each basis function "phi_i" *************************************************/

// Version 1: "u" is an arma::vec 
arma::vec L2project_Vec(	const my_mesh::MeshData &mesh,
				const arma::vec &u)
{
	arma::vec result(mesh.num_nodes);
	result.zeros();

	for (int t=0; t<mesh.num_elements; t++)
	{	std::vector<int> neigh_nodes = mesh.topology2to0[t];
		double area_t = mesh.ele_areas[t];

		for (int i=0; i<neigh_nodes.size(); i++)
		{	int node_i = neigh_nodes[i];
			result(node_i) += area_t / 3.0 * u(node_i);
		}
	}
	return result;
}

// Version 2: "u" is a function
arma::vec L2project_Func(	const my_mesh::MeshData &mesh,
				double (*f)(double, double)	)
{
	arma::vec result(mesh.num_nodes);
	result.zeros();

	for (int t=0; t<mesh.num_elements; t++)
	{	std::vector<int> neigh_nodes = mesh.topology2to0[t];
		double area_t = mesh.ele_areas[t];

		for (int i=0; i<neigh_nodes.size(); i++)
		{	int node_i = neigh_nodes[i];
			double x_i=mesh.nodes[node_i][0], y_i=mesh.nodes[node_i][1];
			result(node_i) += area_t / 3.0 * f(x_i, y_i);
		}
	}
	return result;
}




}
