#include <armadillo>

#include "my_fem.hpp"
#include "../triangle/mesh.hpp"

/* 1.2 C_ij = c* phi_i * phi_j *dx(Omega) */

using namespace arma;
using namespace my_mesh;

namespace my_fem
{

mat assembleMatrixC(const MeshData& mesh, const vec &c)
{
	/* Initialize C to a zero matrix of the right size */
	mat C = zeros(mesh.num_nodes, mesh.num_nodes);

	/* Define matrix entries */
	for (int t=0; t<mesh.num_elements; t++)
	{
		int v0 = mesh.elements[t][0];
		int v1 = mesh.elements[t][1];
		int v2 = mesh.elements[t][2];

		C(v0, v0) += mesh.ele_areas[t] * c(v0) / 3.0;
		C(v1, v1) += mesh.ele_areas[t] * c(v1) / 3.0;
		C(v2, v2) += mesh.ele_areas[t] * c(v2) / 3.0;
	}
	return C;	 
}


}
