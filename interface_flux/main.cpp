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
	std::string filename = "Cosine";
	double y_control = 0.5;

	/* 1. Mesh */
	my_mesh::MeshData mesh;
	mesh = generateMesh_CosInterface(filename, y_control);
	plot_InterfaceSTLVector(mesh, mesh.interface_curvatures, "curvatures");	// plot interface curvature

	// 2. State equation 
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
