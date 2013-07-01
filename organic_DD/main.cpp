#include <iostream>
#include <vector>
#include <string>
#include <armadillo>

#include "../triangle/mesh.hpp"
#include "../gnuplot/my_gnuplot.hpp"
//#include "../my_fem/my_fem.hpp"

#include "gummel_iteration.hpp"
#include "auxillary.hpp"
#include "test.hpp"


int main (int argc, char* argv[])
{
using namespace std;

	// 0. Output some parameters
	std::cout << "\n\n";
	std::cout << "lambda_squared is " << parameters::lambda_squared << "\n";
	std::cout << "epsilon_r = " << parameters::epsilon_rel << "\n";
	std::cout << "mu_n(donor) = " << parameters::mu_n_donor<<"\tmu_n(acceptor) = " << parameters::mu_n_acceptor << "\n";
	std::cout << "mu_p(donor) = " << parameters::mu_p_donor<<"\tmu_p(acceptor) = " << parameters::mu_p_acceptor << "\n";
	std::cout << "mu_x = " << parameters::mu_x << "\n";

	string filename = "organic_DD";

	// 1. Mesh
	vector< vector<double> > interface_nodes;
	vector<double> node(2);
	node[0] = 0.0;	node[1] = 1.0;	interface_nodes.push_back(node);
	node[1] = 0.0;	interface_nodes.push_back(node);

	double max_area = 2e-3;
	my_mesh::MeshData mesh = my_mesh::generateMesh (filename, interface_nodes, max_area);


	// 2. Solve for solution
	// 2.1 solve
	arma::vec psi, n, p, u;
	double applied_psi = -0.0;
	solve_GummelIteration(mesh, psi, n, p, u, applied_psi);

	// 2.2 Plot solution
	bool is_to_plot=true;
	if (is_to_plot)
	{
		my_gnuplot::plot_ArmaVec(mesh, psi, filename+"_psi", "wxt");
//		my_gnuplot::plot_ArmaVec(mesh, n, filename+"_n", "wxt");
//		my_gnuplot::plot_ArmaVec(mesh, p, filename+"_p", "wxt");
//		my_gnuplot::plot_ArmaVec(mesh, u, filename+"_u", "wxt");
		my_gnuplot::plot_ArmaVec(mesh, p-n, filename+"_p-n", "wxt");
	}



	// 3. test
	arma::vec psi_test;
	test_Compute_psi(mesh, n, p, psi_test, applied_psi);
	my_gnuplot::plot_ArmaVec(mesh, psi_test, filename+"_psi_test", "wxt");

	return 0;
}
