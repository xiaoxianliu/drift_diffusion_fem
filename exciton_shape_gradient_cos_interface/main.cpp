#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include <armadillo>

#include "exciton_cos_interface.hpp"
#include "../triangle/mesh.hpp"

#define PI 3.1415926535897932384626433832795

int main()
{
using std::vector;

	std::string filename = "Cosine";


	/********** 1. Mesh */
	// 1.1 Set limit of interface
	double x_control = 1.;			// amplitude of cosine interface
	double x_offset = -0.5;			// x-coordinate of left-most point on interface

	// Determine if "x_control" is valid 
	double x_limit = 0.95;
	if ( x_offset + x_control > x_limit )
	{	std::cout << "x_offset is " << x_offset << "\n";
		std::cout << "x_control (amplitude of cos curve) is " << x_control << "\n";
		std::cout << "abs(x_offset +/- x_control) has to be smaller than " << x_limit << "\n";
		exit(1);
	}

	// 1.2 Generate nodes on interface 
	int num_interface_nodes = 201;
	vector< vector<double> > interface_nodes;

	for (int i=0; i<num_interface_nodes; i++)
	{	vector<double> node(2);
		double y;
		if (i<num_interface_nodes-1)
			y = 1.0 - i/static_cast<double>(num_interface_nodes-1);
		else
			y = 0.0;
		double x = x_offset + x_control * 0.5*(1 - std::cos(2*PI*y));

		node[0] = x;
		node[1] = y;

		interface_nodes.push_back(node);
	}


	// 1.3 Generate mesh
	my_mesh::MeshData mesh;
	double max_area = 0.001;
	mesh = my_mesh::generateMesh(filename, interface_nodes, max_area);
	my_mesh::ComputeInterfaceCurvatures(mesh);				// compute interface curvature
	plot_InterfaceSTLVector(mesh, mesh.interface_curvatures, "curvatures");	// plot interface curvature



	/***** 2. State equation  */
	arma::vec u = solveStateEq(mesh);
	plot_ArmaVec(mesh, u, filename+"_state");
	plot_ArmaVec_on_Interface(mesh, u, filename+"_state");			// plot State Eq solution on interface



	// 3. Adjoint equation 
	arma::vec xi = solveAdjointEq(mesh);
	plot_ArmaVec(mesh, xi, filename+"_adjoint");
	plot_ArmaVec_on_Interface(mesh, xi, filename+"_adjoint");		// plot Adjoint Eq solution on interface

	// 4. Compute shape gradient 
	arma::vec shape_grad = computeShapeGradient(mesh, u, xi);
	plot_InterfaceArmaVec(mesh, shape_grad, "shape_gradient");		// shape gradient



	/* 5. Compute functional */
	double loss = interfaceIntegral( mesh, u );
	std::cout << "total loss at interface is " << loss << "\n";





	return 0;
}
