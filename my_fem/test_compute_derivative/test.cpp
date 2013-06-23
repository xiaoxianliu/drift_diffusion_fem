#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <armadillo>
#include "../../triangle/mesh.hpp"
#include "../my_fem.hpp"
#include "test.hpp"

/* This file computes the gradient of a given function "f" (represented in an FEM space already, of course).



*/



int main(int argc, char* argv[])
{
	std::string filename = "compute_grad";
	// 1. Mesh
	std::vector< std::vector<double> > interface_nodes;
	std::vector<double> new_node(2);
	new_node[0] = 0.0;	new_node[1] = 1.0;	interface_nodes.push_back(new_node);
	new_node[0] = 0.0;	new_node[1] = 0.0;	interface_nodes.push_back(new_node);


	// 2. Compute solution on a successive sequence of meshes
	int level = 3;				// mesh info
	double max_length_0 = 1e-1;
	std::vector<double> lengths;

	std::vector<double> L2_errors;		// error info

	my_mesh::MeshData mesh;
	arma::vec psi_grad_x, psi_grad_y;	// FEM gradient

	for (int i=0; i<level; i++)
	{	double max_length = max_length_0 / pow(2, i);
		lengths.push_back(max_length);

		// 2.1 Generate mesh for current max size
		double max_area = pow(max_length,2);
		mesh = my_mesh::generateMesh (filename, interface_nodes, max_area);

		// 2.2 Interpolate function and its gradient
		arma::vec psi = my_fem::interpolateFunction (mesh, f);
		arma::vec psi_grad_x_interpolated = my_fem::interpolateFunction (mesh, grad_f_x);
		arma::vec psi_grad_y_interpolated = my_fem::interpolateFunction (mesh, grad_f_y);

		// 2.3 Compute gradients by Galerkin method
		my_fem::computeGradient (mesh, psi, psi_grad_x, psi_grad_y);

		// 2.4 Compute error vector
		arma::vec error_x = psi_grad_x - psi_grad_x_interpolated;

		double L2_error;
		L2_error = sqrt (my_fem::integrate_Domain(mesh, error_x % error_x));
		L2_errors.push_back(L2_error);
	}


	// 3. Store L2 error information and plot it

	// 3.1 Dat file
	std::string dat_filename = filename + "_error.dat";
	std::ofstream dat_stream;
	dat_stream.open (dat_filename.c_str());

	for (int i=0; i<level; i++)
	{	dat_stream << lengths[i] << "\t" << L2_errors[i] << "\n";
	}
	dat_stream.close();

	// 3.2 Gnuplot file
	std::string gnuplot_filename = filename + "_error.gnuplot";
	std::ofstream gnuplot_stream;
	gnuplot_stream.open (gnuplot_filename.c_str());
	gnuplot_stream << "plot \"" << dat_filename << "\" with lines\n";
	gnuplot_stream.close();

	// 3.3 Plot
	std::string cmd = "gnuplot \"" + gnuplot_filename + "\" --persist \n";
	system(cmd.c_str());


	// 3.4 Error order
	double order = log(L2_errors[0]/L2_errors[1]) / log(2);
	std::cout << "Order of accuracy is " << order << "\n";


	// 4. Plot solution
	plot_ArmaVec (mesh, psi_grad_x, filename+"_x");
	plot_ArmaVec (mesh, psi_grad_y, filename+"_y");

	return 0;
}



