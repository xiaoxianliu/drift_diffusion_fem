#include <armadillo>
#include "../triangle/mesh.hpp"
#include "misc.hpp"

arma::vec computeShapeGradient(	const my_mesh::MeshData &mesh,
				const arma::vec &u, const arma::vec &xi)
{
	int num_interface_nodes = mesh.interface_nodes.size();

	// 1. The part of gradient involving normal derivative of u
	arma::vec du_dnu_1, du_dnu_2;
	computeInterfaceNormalDerivative(mesh, u, du_dnu_1, du_dnu_2);

	arma::vec g1(num_interface_nodes);
	for (int i=0; i<num_interface_nodes; i++)
	{
		int node_i = mesh.interface_nodes[i];
		g1(i) = (xi(node_i) - 1.0) * ( - du_dnu_2(i) );
	}

	// 2. The part involving normal derivative of xi
	arma::vec dxi_dnu_1, dxi_dnu_2;
	computeInterfaceNormalDerivative(mesh, xi, dxi_dnu_1, dxi_dnu_2);

	arma::vec g2(num_interface_nodes);
	for (int i=0; i<num_interface_nodes; i++)
	{
		int node_i = mesh.interface_nodes[i];
		g2(i) = dxi_dnu_1(i) * u(node_i);
	}

	// 3. The part involving curvature
	arma::vec g3(num_interface_nodes);
	for (int i=0; i<num_interface_nodes; i++)
	{
		int node_i = mesh.interface_nodes[i];
		g3(i) = (xi(node_i) - 1.0) * u(node_i) * mesh.interface_curvatures[i];
	}

	arma::vec g;
	g = g1 + g2 + g3;

	return g;
}
