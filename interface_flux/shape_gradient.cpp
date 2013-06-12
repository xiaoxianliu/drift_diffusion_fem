#include <vector>

#include <armadillo>

#include "../triangle/mesh.hpp"

#include "misc.hpp"			// include functions related to the definition of state and adjoint equations




// Declare function to be used in "computeShapeGradient" 
double gradv_dot_gradu_over_one_triangle(	const my_mesh::MeshData &mesh, 
						const int t,		// "t" being the index of the chosen element in "mesh"
						const arma::vec &u,
						const arma::vec &v);




/*******************************************************************************/
/***** compute shape gradient (function to be included in main.cpp) ************/

arma::vec computeShapeGradient(	const my_mesh::MeshData &mesh,
				const arma::vec u,					// solution to state equation
				const arma::vec xi)					// solution to adjoint equation
{
	int num_interface_nodes = mesh.interface_nodes.size();
	int num_interface_edges = num_interface_nodes - 1;






	// M * g = rhs

	// 1. Define coefficient matrix, which is diagonal because of 1st-order trapozoidal-rule approximation

	arma::mat M(num_interface_nodes, num_interface_nodes);	// square coefficient matrix of dimension "num_interface_nodes"
	M.zeros();						// initialized with zero matrix

	for (int i=0; i<num_interface_edges; i++)	// recall mesh.interface_nodes and mesh.interface_edges are in stored in order
	{
		int edge_global_index = mesh.interface_edges[i];
		double edge_length = mesh.edge_lengths[edge_global_index];

		M(i,i) += edge_length / 2.0;
		M(i+1, i+1) += edge_length / 2.0;
	}





	// 2. rhs involves 3 parts
	// 2.1 the part involving value on interface, including curvature
	arma::vec rhs0(num_interface_nodes);
	rhs0.zeros();

	for (int i=0; i<num_interface_nodes; i++)
	{	int node_global_index = mesh.interface_nodes[i];
		rhs0(i) = ( xi( node_global_index ) - 1 ) * u( node_global_index ) * mesh.interface_curvatures[i];
	}



















	// 2.2 The part involving normal derivative of region "1"

	arma::vec rhs1(num_interface_nodes);
	rhs1.zeros();

	for (int i=0; i<num_interface_nodes; i++)
	{
		int node_global_index = mesh.interface_nodes[i];


		arma::vec phi_u(mesh.num_nodes);					// arma::vec representation of "phi * u"
		phi_u.zeros();
		phi_u(node_global_index) = u(node_global_index);


		std::vector<int> neigh_elements = mesh.topology0to2[ node_global_index ];	// indices of neighboring elements
												// to "node_global_index"-th node
//		std::cout << "this is the " << i << "-th node on interface with global index " << node_global_index << "\n";

		for (int j=0; j<neigh_elements.size(); j++)
		{	int t = neigh_elements[j];				// global_index for each neighboring triangle

//			std::cout << j << "-th neighboring element has global index " << t << "\n";		

			if (mesh.element_markers[t] == 1)
			{
				rhs1(i) += mesh.ele_areas[t] / 3.0 * u(node_global_index) * xi(node_global_index);
				rhs1(i) += gradv_dot_gradu_over_one_triangle(mesh, t, phi_u, xi);
			}
		}


	}








	// 2.3 the part involving normal derivative of region "2"
	arma::vec rhs2(num_interface_nodes);
	rhs2.zeros();

	for (int i=0; i<num_interface_nodes; i++)
	{
		int node_global_index = mesh.interface_nodes[i];
		double x = mesh.nodes[node_global_index][0];
		double y = mesh.nodes[node_global_index][1];

		arma::vec xi_phi(mesh.num_nodes);			// xi_phi denote the linear interpolation of "(xi-1)*phi"
		xi_phi.zeros();
		xi_phi(node_global_index) = xi(node_global_index) - 1 ;

		std::vector<int> neigh_elements = mesh.topology0to2[ node_global_index ];	// neighboring elements to this node
		for (int j=0; j<neigh_elements.size(); j++)
		{
			int t = neigh_elements[j];

			if ( mesh.element_markers[t] == 2 )
			{
				rhs2(i) += mesh.ele_areas[t]/3.0 
					  * ( xi(node_global_index) - 1.0 )
					  * ( func_g(x,y) - u(node_global_index));
				rhs2(i) -= gradv_dot_gradu_over_one_triangle(mesh, t, xi_phi, u);
			}
		}
	}


	// Finally add all 3 "rhs" together
	arma::vec rhs = rhs0 
			+ rhs1
			+ rhs2
			;









	// Now solve for shape gradient function

	arma::vec g(num_interface_nodes);					// initialize shape gradient vector "g"
	g.zeros();
	g = arma::solve(M, rhs);


	return g;
}





















/* Start of definition of function to be called by "computeShapeGradient(...)" */


/* "u" and "v" are two linear function defined over "mesh". */
/* "gradv_dot_gradu"(...) computes the dot product of grad(v) and grad(u) over a given triangular element with index "ele_index" */
double gradv_dot_gradu_over_one_triangle(	const my_mesh::MeshData &mesh, 
						const int t,		// "t" being the index of the chosen element in "mesh"
						const arma::vec &u,
						const arma::vec &v)
{
	int v0 = mesh.elements[t][0];
	int v1 = mesh.elements[t][1];
	int v2 = mesh.elements[t][2];

	/* Compute J_inv */
	arma::vec r1(2);						// difference vector: node_1 - node_0
	arma::vec r2(2);						// difference vector: node_2 - node_0
	r1(0) = mesh.nodes[v1][0] - mesh.nodes[v0][0];	r1(1) = mesh.nodes[v1][1] - mesh.nodes[v0][1];
	r2(0) = mesh.nodes[v2][0] - mesh.nodes[v0][0];	r2(1) = mesh.nodes[v2][1] - mesh.nodes[v0][1];

	arma::mat J(2,2);
	J(0,0) = arma::dot(r1, r1);	J(0,1) = arma::dot(r1, r2);
	J(1,0) = arma::dot(r1, r2);	J(1,1) = arma::dot(r2, r2);
	arma::mat J_inv = arma::inv(J);

	/* arma::vec for nodal differences in both "u" and "v" */
	arma::vec delta_u(2);
	arma::vec delta_v(2);
	delta_u(0) = u(v1) - u(v0);	delta_u(1) = u(v2) - u(v0);
	delta_v(0) = v(v1) - v(v0);	delta_v(1) = v(v2) - v(v0);

	/* Compute <grad_u, grad_v> */
	double result;
	result = arma::dot( J_inv*delta_u, delta_v);

	return result;
}

/* End of definition of function to be called by "computeShapeGradient(...)" */

