#include <armadillo>
#include "../triangle/mesh.hpp"

using namespace arma;
using namespace my_mesh;

namespace my_fem
{

mat assembleMatrixD(const MeshData& mesh, const vec &d)
{
	mat D(mesh.num_nodes, mesh.num_nodes);
	D.zeros();

	for (int e=0; e<mesh.num_edges; e++)
	{	if (mesh.edge_markers[e]==5)
		{
			int v0 = mesh.edges[e][0];
			int v1 = mesh.edges[e][1];
			double length = mesh.edge_lengths[e];

			D(v0, v0) += d(v0) * length / 2.0;
			D(v1, v1) += d(v1) * length / 2.0;
		}
	}
	return D;
}

}

