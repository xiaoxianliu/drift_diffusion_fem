#include <iostream>
#include <cstdlib>
#include <vector>
#include <armadillo>

#include "../triangle/mesh.hpp"
#include "misc.hpp"


// Declare local function(s)
double dphi_dx_dot_du_dx(const my_mesh::MeshData &mesh, const arma::vec &u, int node_i, int ele_i);


/******** Main function in this file *************/
int computeInterfaceNormalDerivative(	const my_mesh::MeshData &mesh, 
					const arma::vec &u,
					arma::vec &du_dnu_1, arma::vec &du_dnu_2) //"nu" indicates the outer normal direction on surface
{
	// Reset "du_dnu_1" and "du_dnu_2" to zero vectors in case a mis-initialization was in place
	int num_interface_nodes = mesh.interface_nodes.size();
	du_dnu_1.zeros( num_interface_nodes );
	du_dnu_2.zeros( num_interface_nodes );


	// assemble linear system by Galerkin method on interface nodes;
	// Left side has regional marker "1"; right side has regional marker "2"

	// 1. Coefficient matrix
	arma::mat M1(num_interface_nodes, num_interface_nodes), M2(num_interface_nodes, num_interface_nodes);
	M1.zeros();
	M2.zeros();

	for (int i=0; i<num_interface_nodes-1; i++)
	{
		double length = mesh.edge_lengths[ mesh.interface_edges[i] ];
		M1(i,i) += 0.5 * length;	M1(i+1, i+1) += 0.5 * length;
		M2(i,i) += 0.5 * length;	M2(i+1, i+1) += 0.5 * length;
	}


	// 2. right-hand side vector
	arma::vec rhs1(num_interface_nodes), rhs2(num_interface_nodes);
	rhs1.zeros();
	rhs2.zeros();

	for (int i=0; i<num_interface_nodes; i++)
	{
		int node_i = mesh.interface_nodes[i];
		double x,y;
		x = mesh.nodes[node_i][0];	y = mesh.nodes[node_i][1];

		// "node_i"-th entry of "rhs" vector
		std::vector<int> neigh_elements = mesh.topology0to2 [node_i];
		for (int dummy=0; dummy < neigh_elements.size(); dummy++)
		{
			int t = mesh.topology0to2[node_i][dummy];

			if ( mesh.element_markers[t] == 1)
			{	double area_t = mesh.ele_areas[t];

				rhs1(i) += area_t * dphi_dx_dot_du_dx(mesh, u, node_i, t);
				rhs1(i) += area_t * ( u(node_i) - func_g(x,y) ) /3.0;
			}
			else if ( mesh.element_markers[t] == 2)
			{	double area_t = mesh.ele_areas[t];

				rhs2(i) += area_t * dphi_dx_dot_du_dx(mesh, u, node_i, t);
				rhs2(i) += area_t * ( u(node_i) - func_g(x,y) ) /3.0;
			}
		}
	}


	du_dnu_1 = arma::solve(M1, rhs1);
	du_dnu_2 = arma::solve(M2, rhs2);

	return 0;
}

















double dphi_dx_dot_du_dx(	const my_mesh::MeshData &mesh,
				const arma::vec &u,
				int node_i, int ele_i)		// "node_i" is the index of the node correponding to basis "phi_i"
								// "ele_i" is the index of the element containing "node_i"
{
	std::vector<int> node_indices = mesh.topology2to0[ele_i];
	int node_j=-1, node_k = -1;

	// Identify the other "2" nodes of element "ele_i"
	for (int i=0; i<3; i++)
	{
		if ( node_indices[i] != node_i && node_j==-1)
			node_j = node_indices[i];
		else if ( node_indices[i] != node_i && node_k==-1)
			node_k = node_indices[i];
	}
	if (node_j==-1 || node_k==-1)
	{	std::cout << "Node_i is " << node_i << ", node_j is " << node_j << ", node_k is " << node_k << "\n";
		std::cout << "All of the 3 should be non-negative. Exit...\n";
		exit(1);
	}

	// Jacobian matrix of transformation from "ele_i" to reference triangle
	arma::mat jacobian(2,2);
	jacobian(0,0) = mesh.nodes[node_j][0] - mesh.nodes[node_i][0];
	jacobian(1,0) = mesh.nodes[node_j][1] - mesh.nodes[node_i][1];
	jacobian(0,1) = mesh.nodes[node_k][0] - mesh.nodes[node_i][0];
	jacobian(1,1) = mesh.nodes[node_k][1] - mesh.nodes[node_i][1];

	// Compute matrix " J^(-1)*J^(-T) "
	arma::mat JT_J(2,2);
	JT_J = jacobian.t() * jacobian;
	arma::mat M = JT_J.i();

	// Finally we're ready to compute < dphi/dx, du/dx > on element "ele_i":
	//	< dphi/dx, du/dx > = < dphi/d_xi, M * du/d_xi >
	arma::vec dphi_d_xi(2);
	dphi_d_xi.ones();	dphi_d_xi *= -1.0;
	arma::vec du_d_xi(2);
	du_d_xi(0) = u(node_j) - u(node_i);	du_d_xi(1) = u(node_k) - u(node_i);

	double result;
	result = arma::dot( dphi_d_xi, M * du_d_xi);

	return result;

}
