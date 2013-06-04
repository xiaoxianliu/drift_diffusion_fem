#include <iostream>
#include <string>

#include <armadillo>

#include "exciton_single_geometry.hpp"
#include "../fem_assemble/fem_assemble.hpp"

int main()
{
	std::string filename = "smart";

	/* 1. Generate mesh */
	my_mesh::MeshData mesh = generateMesh("smart");

	/* 2. Form coefficient matrix and RHS vector */
	/* 2.1 Coefficient matrix */
	arma::vec a = interpolate_func(mesh, diffusion_coefficient);
	arma::vec c = interpolate_func(mesh, decay_rate);
	arma::vec d = interpolate_func(mesh, interface_reaction_rate);

	arma::mat A = linear_fem::assembleMatrixA(mesh, a);
	arma::mat C = linear_fem::assembleMatrixC(mesh, c);
	arma::mat D = linear_fem::assembleMatrixD(mesh, d);

	arma::mat M = A + C + D;

	/* 2.2 RHS vector */
	arma::vec g = project_func(mesh, generation_rate);

	std::cout<<"g is \n" << g << "\n";


	return 0;
}
