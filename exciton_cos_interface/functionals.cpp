#include <iostream>
#include <armadillo>
#include "../triangle/mesh.hpp"

#include "misc.hpp"


/* y = int_Gamma f*u*ds */
double interfaceIntegral(const my_mesh::MeshData& mesh, const arma::vec& u )
{	double sum=0;


	for (int i=0; i < mesh.interface_edges.size(); i++)
	{	int node0 = mesh.interface_nodes[i];			// two end nodes of interface
		int node1 = mesh.interface_nodes[i+1];

		double x0 = mesh.nodes[node0][0];
		double y0 = mesh.nodes[node0][1];
		double x1 = mesh.nodes[node1][0];
		double y1 = mesh.nodes[node1][1];

		double length = mesh.edge_lengths[ mesh.interface_edges[i] ];

		double increment = 0.5 * length * ( func_d(x0,y0)*u(node0) + func_d(x1,y1)*u(node1) );

		std::cout << "value of u is " << u(node0) << "\n";
		std::cout << "length is " << length << "\n";
		std::cout << "increment is " << increment << "\n";


		sum += increment;

	}

/*	for (int e=0; e<mesh.num_edges; e++)
	{
		if (mesh.edge_markers[e]==5)
		{	int v0 = mesh.edges[2*e];
			int v1 = mesh.edges[2*e+1];
			double length = mesh.edge_lengths[e];

			sum += length * (u(v0)*d(v0) + u(v1)*d(v1)) / 2.0;		// Trapezoidal rule for 1st order approx.
		}
	}
*/
	return sum;
}
