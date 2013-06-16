#include <armadillo>
#include "../triangle/mesh.hpp"
#include "fem_assemble.hpp"

using namespace my_mesh;
using namespace arma;

namespace linear_fem
{

mat assembleMatrixA(const MeshData &mesh, const vec &a)
{
	int num_nodes = mesh.num_nodes;

	/* Initialize matrix to all 0's*/
	mat A=arma::zeros(num_nodes, num_nodes);

	/* Define matrix entries */
	for (int t=0; t<mesh.num_elements; t++)			// loop through all elments
	{	if (mesh.num_nodes_per_ele!=3)			// 3 for linear element; 6 for quadratic element, etc
			{std::cout << "Element has to be linear... i.e. 3 nodes per element " << std::endl; exit(1);}
		/* Vertices of triangle "t" */
		int v0 = mesh.elements[t][0];
		int v1 = mesh.elements[t][1];
		int v2 = mesh.elements[t][2];

		vec r0(2), r1(2), r2(2);			// 2-dim vectors of coordinates;
		r0(0) = mesh.nodes[v0][0];	r0(1) = mesh.nodes[v0][1];
		r1(0) = mesh.nodes[v1][0];	r1(1) = mesh.nodes[v1][1];
		r2(0) = mesh.nodes[v2][0];	r2(1) = mesh.nodes[v2][1];

		/* factor shared by all "A_ij(t)" */
		double factor;
		factor = (a(v0) + a(v1) + a(v2))/3.0 /(4.0*mesh.ele_areas[t]);

		/* Modify the 9 entries corresponding to the 3 vertices */
		A(v0,v0) += factor * dot(r1-r2, r1-r2);
		A(v0,v1) += factor * (- dot(r0-r2, r1-r2));
		A(v0,v2) += factor * (- dot(r2-r1, r0-r1));

		A(v1,v1) += factor * dot(r0-r2, r0-r2);
		A(v1,v0) += factor * (- dot(r0-r2, r1-r2));
		A(v1,v2) += factor * (- dot(r1-r0, r2-r0));

		A(v2,v2) += factor * dot(r0-r1, r0-r1);
		A(v2,v1) += factor * (- dot(r1-r0, r2-r0));
		A(v2,v0) += factor * (- dot(r2-r1, r0-r1));
	}
	return A;
}




}
