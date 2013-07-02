#include <iostream>
#include <cmath>
#include <vector>

#include <armadillo>

#include "../triangle/mesh.hpp"
#include "../my_fem/my_fem.hpp"
#include "parameters.hpp"

#define PI 3.14159265359

//// Auxillary functions that is needed by equation solvers

/*******************************************************************************************************************/
// Compute magnitude of electrical field |E|
int compute_ElectricFieldAmplitude(	const my_mesh::MeshData &mesh,
					const arma::vec &psi,
					arma::vec &E)
{
	// 1. compute electrical field vector by components
	arma::vec grad_psi_x, grad_psi_y;
	my_fem::computeGradient(mesh, psi, grad_psi_x, grad_psi_y);

	// 2. compute the magnitude (i.e. non-negative) of electrical field
	E.resize(mesh.num_nodes);
	for (int i=0; i<mesh.num_nodes; i++)
	{
		E(i) = sqrt( pow(grad_psi_x(i), 2) + pow(grad_psi_y(i), 2) );
	}
	return 0;
}















/*********************************************************************************************************************/
// mobility of electrons and holes
// 1. element-wise mobility
// 1.1 electron
int	compute_MobilityN_elementwise(	const my_mesh::MeshData &mesh,
					const arma::vec &E,				// electric field intensity
					arma::vec &mu_n_elementwise)
{
	// 1. Copy the corresponding parameters from namespace "parameters"
	double mu_n_1 = parameters::mu_n_donor;
	double mu_n_2 = parameters::mu_n_acceptor;
	double gamma = parameters::gamma_n;

	// 2. Compute average mobility for each element
	for (int t=0; t<mesh.num_elements; t++)
	{	int v0 = mesh.topology2to0[t][0];
		int v1 = mesh.topology2to0[t][1];
		int v2 = mesh.topology2to0[t][2];

		// Determine zero-field mobility of electron
		double mu_n_t_zerofield;
		if (mesh.element_markers[t]==1)
			mu_n_t_zerofield = mu_n_1;
		else if (mesh.element_markers[t] == 2)
			mu_n_t_zerofield = mu_n_2;
		else
		{	std::cout << "marker of element " << t << " is " << mesh.element_markers[t] << ". It has to be 1 or 2.\n";
			exit(1);
		}

		mu_n_elementwise(t) = mu_n_t_zerofield *
					  ( exp(gamma * sqrt(E(v0))) 
					   +exp(gamma * sqrt(E(v1)))
					   +exp(gamma * sqrt(E(v2)))
					  )/3.0;
	}

	return 0;
}

// 1.2 hole
int compute_MobilityP_elementwise (	const my_mesh::MeshData &mesh, 
					const arma::vec &E, 
					arma::vec &mu_p_elementwise)
{
	// 1. Copy the corresponding parameters from namespace "parameters"
	double mu_p_1 = parameters::mu_p_donor;
	double mu_p_2 = parameters::mu_p_acceptor;
	double gamma_p = parameters::gamma_p;

	// 2. Compute average mobility for each element
	for (int t=0; t<mesh.num_elements; t++)
	{	int v0 = mesh.topology2to0[t][0];
		int v1 = mesh.topology2to0[t][1];
		int v2 = mesh.topology2to0[t][2];

		// Determine zero-field mobility of electron
		double mu_p_t_zerofield;
		if (mesh.element_markers[t]==1)
			mu_p_t_zerofield = mu_p_1;
		else if (mesh.element_markers[t] == 2)
			mu_p_t_zerofield = mu_p_2;
		else
		{	std::cout << "marker of element " << t << " is " << mesh.element_markers[t] << ". It has to be 1 or 2.\n";
			exit(1);
		}

		mu_p_elementwise(t) = mu_p_t_zerofield *
					  ( exp(gamma_p * sqrt(E(v0))) 
					   +exp(gamma_p * sqrt(E(v1)))
					   +exp(gamma_p * sqrt(E(v2)))
					  )/3.0;
	}

	return 0;
}

// 1.3 Exciton mobility
int compute_MobilityX_elementwise (	const my_mesh::MeshData &mesh,
					arma::vec &mu_x_elementwise)
{
	double mu_x_1 = parameters::mu_x_donor;
	double mu_x_2 = parameters::mu_x_acceptor;

	for (int t=0; t<mesh.num_elements; t++)
	{	if (mesh.element_markers[t]==1)
			mu_x_elementwise(t) = mu_x_1;
		else if (mesh.element_markers[t] == 2)
			mu_x_elementwise(t) = mu_x_2;
	}
return 0;
}













/*******************************************************************************************************************/
/* Compute "outward" normal derivative of given vector "u" on interface w.r.t. either left region (marked "1") or right region (marked "2") */

// Local function used by "compute_dpsi_dnu_interface(...)"
// Input: linear interpolated function "u", linear basis function at node_i "phi_i", chosen element index "ele_i"
// Output: grad(phi_i)_dot_grad(u) on triangle element "ele_i"; it's a piecewise constant on different elements
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






/////// Main function
int compute_dpsi_dnu_interface(	const my_mesh::MeshData &mesh,
					const arma::vec &psi,
					const arma::vec &n,
					const arma::vec &p,
					const arma::vec &u,
					arma::vec &dpsi_dnu_1,
					arma::vec &dpsi_dnu_2)
{
	int num_interface_nodes = mesh.interface_nodes.size();
	double epsilon1 = parameters::epsilon_donor;
	double epsilon2 = parameters::epsilon_acceptor;
	double lambda_squared = parameters::lambda_squared;

	// 1. Coefficient
	arma::mat M(num_interface_nodes, num_interface_nodes);
	M.zeros();

	for (int i=0; i < num_interface_nodes-1 ; i++)		// recall num_interface_edges = num_interface_nodes - 1
	{	int v0 = mesh.interface_nodes[i];
		int v1 = mesh.interface_nodes[i+1];
		int edge_index = mesh.interface_edges[i];
		double length = mesh.edge_lengths [edge_index];
		M(i,i) += 0.5 * length;
		M(i+1, i+1) += 0.5*length;
	}

	// 2. right-hand-side vectors
	arma::vec rhs1(num_interface_nodes), rhs2(num_interface_nodes);
	rhs1.zeros();	rhs2.zeros();
	for (int i=0; i < num_interface_nodes; i++)
	{	int node_index = mesh.interface_nodes[i];
		std::vector<int> neigh_elements = mesh.topology0to2[node_index];	// neighboring elements of node[node_index]

		for (int t=0; t<neigh_elements.size(); t++)
		{	int element_index = neigh_elements[t];
			double area_t = mesh.ele_areas[element_index];

		// for each triangle "t", compute "int_t [grad^2(psi)*phi_i + grad(psi)*grad(phi_i)]*dx
			if ( mesh.element_markers [element_index] == 1)
			{	double increment0, increment1;
				increment0 = ( (n(node_index) - p(node_index))/(lambda_squared * epsilon1) )*area_t/3.0;
				increment1 = dphi_dx_dot_du_dx (mesh, psi, node_index, element_index) *area_t;
				rhs1(i) += increment0 + increment1;
			}
			if ( mesh.element_markers [element_index] == 2)
			{	double increment0, increment1;
				increment0 = ( (n(node_index) - p(node_index))/(lambda_squared * epsilon2) )*area_t/3.0;
				increment1 = dphi_dx_dot_du_dx (mesh, psi, node_index, element_index) *area_t;
				rhs2(i) += increment0 + increment1;
			}
		}
	}

	dpsi_dnu_1 = arma::solve(M, rhs1);
	dpsi_dnu_2 = arma::solve(M, rhs2);
	
	return 0;
}





/*******************************************************************************************************************/
/* Compute exciton dissociation rate.
 Only nonzero on the interface, but defined on all nodes */

// Local function compute enhancement factor of dissociation rate for given "dpsi_dnu_1" on some node
// Ref: 2003 Barker, Ramsdale, and Greenham PRB. 67. 075205
double electric_field_factor_Func(double dpsi_dnu_1)
{
using namespace parameters;
	double result;

	double M=0;
	M = sqrt( q_unit * fabs(dpsi_dnu_1)/(PI*epsilon_vacuum * epsilon_rel) )/ U_T; 
	if ( M<1e-6 )
		result = 1.0;				// almost zero electrical field
	else if (dpsi_dnu_1>=0)	
		result = 2.0 * (exp(M)*(1.-1./M) + 1/M) / M;		// "dpsi_dnu_1 >0" is equivalent to "E_nu_1 <0"
	else
		result = 4.0 / pow(M,2) * ( 1 - exp(-pow(M,2)/4) );
	return result;

}


//////  Main function
int compute_ExcitonDissociationRate(	const my_mesh::MeshData &mesh,
					const arma::vec &psi,
					const arma::vec &n,
					const arma::vec &p,
					const arma::vec &u,
					arma::vec &k_diss_interface)
{
	// 1. Compute dpsi/dnu on interface
	arma::vec dpsi_dnu_1, dpsi_dnu_2;
	compute_dpsi_dnu_interface (mesh, psi, n, p, u, dpsi_dnu_1, dpsi_dnu_2);

	// 2. Compute dissociation rate node by node; it's only nonzero on the interface
	k_diss_interface.resize(mesh.num_nodes);
	k_diss_interface.zeros();

	using parameters::h_interface;						// interface width
	using parameters::k_diss_0;						// zero-field dissociation rate

	for (int i=0; i<mesh.interface_nodes.size(); i++)
	{	int node_index = mesh.interface_nodes[i];
		k_diss_interface(node_index) = h_interface * k_diss_0 * electric_field_factor_Func( dpsi_dnu_1(i) );
	}

	return 0;
}














/*********************************************************************************************************************/
// Recombination rate
int compute_RecombinationRate(	const my_mesh::MeshData &mesh,
				const arma::vec &E,				// electric field amplitude
				arma::vec &recomb_interface)
{
	// 1. Reset recomb_interface to the correct dimension with 0's
	recomb_interface.resize(mesh.num_nodes);
	recomb_interface.zeros();

	// 2. Identify the mobility to use on interface
	double mu_n_0 = parameters::mu_n_donor > parameters::mu_n_acceptor ? parameters::mu_n_donor : parameters::mu_n_acceptor;
	double mu_p_0 = parameters::mu_p_donor > parameters::mu_p_acceptor ? parameters::mu_p_donor : parameters::mu_p_acceptor;

	// 3. Compute recombination rate on interface
	using parameters::h_interface;						// dimensionless interface width
	using parameters::gamma_n;
	using parameters::gamma_p;
	using parameters::recomb_coefficient;
	using parameters::epsilon_rel;

	for (int i=0; i<mesh.interface_nodes.size(); i++)
	{	int node_index = mesh.interface_nodes[i];
		double volume_recomb = recomb_coefficient * (	  mu_n_0 * exp(gamma_n * sqrt(E(node_index)))
								+ mu_p_0 * exp(gamma_p * sqrt(E(node_index)))
							      ) / epsilon_rel;	// volume recombination rate
		double interface_recomb = volume_recomb * h_interface;		// surface recombination rate
		recomb_interface(node_index) = interface_recomb;
	}
	return 0;
}








/*************************************************************************************************************************/
// Compute arma::vec of "photo generation function"
double photo_generation_Func(double x, double y)
{
using namespace parameters;
	double x_incident = -1.0;
	return Q0*exp(-alpha*(x - x_incident));		// assuming light comes from the left
}

int compute_PhotoGenerationVec(	const my_mesh::MeshData &mesh,
				arma::vec &Q_vec)
{
	Q_vec = my_fem::interpolateFunction(mesh, photo_generation_Func);
return 0;
}





