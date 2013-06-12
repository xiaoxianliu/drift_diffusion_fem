#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include <armadillo>

#include "exciton_interface_flux.hpp"
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


	// 3. Compute normal derivative on interface
	arma::vec du_dnu_1, du_dnu_2;
	computeInterfaceNormalDerivative (mesh, u, du_dnu_1, du_dnu_2);
	plot_InterfaceArmaVec(mesh, du_dnu_1, filename+"_interface_normal_deriv_1");
	plot_InterfaceArmaVec(mesh, du_dnu_2, filename+"_interface_normal_deriv_2");



	/* 5. Compute functional */
	double loss = interfaceIntegral( mesh, u );
	std::cout << "total loss at interface is " << loss << "\n";





	return 0;
}
