#include <vector>
#include <armadillo>
#include "../triangle/mesh.hpp"
#include "exciton_smart_geometry.hpp"


int applyDirichletBC(const my_mesh::MeshData& mesh, arma::mat& M, arma::vec& b)
{

using namespace std;
using namespace arma;

	int num_nodes = mesh.num_nodes;
	for (int v=0; v<num_nodes; v++)
	{
		vector<int> Es = mesh.topology0to1[v];
		for (int i=0; i<Es.size(); i++)
		{	int edge = Es[i];
			int edge_marker = mesh.edge_markers[edge];

			if (edge_marker == 1)			// left boundary condition
			{	double x = mesh.nodes[v][0];
				double y = mesh.nodes[v][1];

				M(v, span::all) = zeros<mat>(1, num_nodes);
				M(v,v) = 1.0;
				b(v) = uD_left(x,y);
			}
			else if (edge_marker == 3)		// right boundary condition
			{	double x = mesh.nodes[v][0];
				double y = mesh.nodes[v][1];

				M(v, span::all) = zeros<mat>(1, num_nodes);
				M(v,v) = 1.0;
				b(v) = uD_right(x,y);
			}
		}

	}

	return 0;
}
