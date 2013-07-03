#include <iostream>
#include <vector>

#include <armadillo>

#include "../../triangle/mesh.hpp"
#include "../../my_fem/my_fem.hpp"
#include "../my_gnuplot.hpp"

double f(double x, double y) { return 1+x+y;}

int main(int argc, char* argv[]){

	// 1. Build mesh
	std::string filename = "test_plot";

	std::vector< std::vector<double> > interface_nodes;
	std::vector<double> new_node(2);
	new_node[0] = 0.0;	new_node[1]= 0.0;	interface_nodes.push_back(new_node);
	new_node[0] = 0.0;	new_node[1]= 1.0;	interface_nodes.push_back(new_node);

	double max_area = 1e-2;

	my_mesh::MeshData mesh;
	mesh = my_mesh::generateMesh(filename, interface_nodes, max_area);

	// 2. Interpolate a given function
	arma::vec f_vec = my_fem::interpolateFunction(mesh, f);

	// 3. Plot functions
	my_gnuplot::plot_ArmaVec(mesh, f_vec, filename);
	my_gnuplot::plot_ArmaVec_on_DirichletBoundary(mesh, 1, f_vec, filename);
	my_gnuplot::plot_ArmaVec_on_DirichletBoundary(mesh, 3, f_vec, filename);
	my_gnuplot::plot_ArmaVec_on_DirichletBoundary(mesh, 2, f_vec, filename);

 return 0;
}
