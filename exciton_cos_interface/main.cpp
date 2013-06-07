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
	std::string filename = "cos_interface";
	double y_control = 0.5;

	/* 1. Mesh */
	my_mesh::MeshData mesh;
	mesh = generateMesh_CosInterface(filename, y_control);

	/* 2. State equation */
	arma::vec u = solveStateEq(mesh);

	/* 3. Adjoint equation */
	arma::vec xi = solveAdjointEq(mesh);

	/* 3. Plot solution to state equaiton */
	plot_ArmaVec(mesh, u, filename+"_state");
	plot_ArmaVec(mesh, xi, filename+"_adjoint");

	plot_ArmaVec_Interface(mesh, u, filename+"_state");			// plot State Eq solution on interface
	plot_ArmaVec_Interface(mesh, xi, filename+"_adjoint");			// plot Adjoint Eq solution on interface

	plot_STLVector_Interface(mesh, mesh.interface_curvatures, "curvatures");// plot interface curvature

	return 0;
}
