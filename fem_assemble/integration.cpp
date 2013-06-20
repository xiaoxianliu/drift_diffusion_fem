#include <vector>
#include <armadillo>
#include "../triangle/mesh.hpp"
#include "fem_assemble.hpp"

namespace linear_fem
{
double integrate_Domain(const my_mesh::MeshData &mesh,			//integration of a piecewise linear function over a 2D domain 
			const arma::vec u)
{
	double sum=0;

	for (int t=0; t<mesh.num_elements; t++)
	{	double area_t = mesh.ele_areas[t];
		std::vector<int> neigh_nodes = mesh.topology2to0[t];
		double nodal_value;					// store nodal value of vector "u", i.e. u(node_i)

		for (int i=0; i<neigh_nodes.size(); i++)		// 3 nodes associated with triangle "t"
		{	int node_i = neigh_nodes[i];
			sum += area_t * u(node_i) / 3.0;
		}
	}
	return sum;
}

}
