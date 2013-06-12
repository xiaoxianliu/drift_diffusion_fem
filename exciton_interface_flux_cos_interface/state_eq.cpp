#include <vector>
#include <iostream>

#include <armadillo>
#include "../triangle/mesh.hpp"
#include "../fem_assemble/fem_assemble.hpp"

#include "exciton_interface_flux.hpp"

#include "misc.hpp"			// include functions related to the definition of state and adjoint equations


/* Main solver for state equation */
arma::vec solveStateEq(const my_mesh::MeshData mesh)
{
using namespace arma;
using std::cout;

	/* 1. Assemble linear system (without imposing Dirichlet BC)*/
	/*	M*u = vec_rhs	*/

	/* 1.1 Coefficient matrix */
	arma::mat M;
	{	arma::vec vec_a = linear_fem::interpolateFunction(mesh, func_a);
		arma::mat A = linear_fem::assembleMatrixA(mesh, vec_a);
//		cout << "vec_a is\n" << vec_a << "\n";
//		cout << "matrix A is \n" << A << "\n";

		arma::vec vec_c = linear_fem::interpolateFunction(mesh, func_c);
		arma::mat C = linear_fem::assembleMatrixC(mesh, vec_c);
//		cout << "vec_c is\n" << vec_c << "\n";
//		cout << "matrix C is \n" << C << "\n";

		arma::vec vec_d = linear_fem::interpolateFunction(mesh, func_d);
		arma::mat D = linear_fem::assembleMatrixD(mesh, vec_d);
//		cout << "vec_d is\n" << vec_d << "\n";
//		cout << "matrix D is \n" << D << "\n";

		M = A + C + D;
	}

	/* 1.2 right-hand side vector */
	arma::vec vec_rhs(mesh.num_nodes);
	vec_rhs.zeros();

	for (int i=0; i<mesh.num_nodes; i++)
	{
		double x = mesh.nodes[i][0];
		double y = mesh.nodes[i][1];
		double g_i = func_g(x,y);

		std::vector<int> neigh_elements = mesh.topology0to2[i];
		for (int j=0; j < neigh_elements.size(); j++)
		{	double area = mesh.ele_areas[ neigh_elements[j] ];
			vec_rhs(i) += area * g_i / 3.0;
		}
	}

//	cout << "vec_rhs is:\n" << vec_rhs << "\n";

	/* 2. Apply Dirichlet boundary condition */
	int num_nodes = mesh.num_nodes;
	for (int i=0; i<num_nodes; i++)
	{
		std::vector<int> Es = mesh.topology0to1[i];			// all neighboring edges to i-th vertex
		for (int j=0; j<Es.size(); j++)
		{	int edge_index = Es[j];
			if ( mesh.edge_markers[edge_index] == 1 || mesh.edge_markers[edge_index] == 3 )
			{	double x = mesh.nodes[i][0];
				double y = mesh.nodes[i][1];
				M(i, arma::span::all) = arma::zeros<arma::mat>(1, num_nodes);
				M(i,i) = 1.0;
				vec_rhs(i) = uD(x,y);

				break;
			}
		}
	}

 	
	/* 3. Solve for solution */
	arma::vec u;
	u = arma::solve(M, vec_rhs);

	return u;
}


