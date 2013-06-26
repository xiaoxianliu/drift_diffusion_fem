#include <iostream>
#include <vector>
#include <string>
#include <armadillo>

#include "../triangle/mesh.hpp"
#include "../gnuplot/my_gnuplot.hpp"

//#include "parameters.hpp"
#include "main.hpp"


int main (int argc, char* argv[])
{
using namespace std;

	string filename = "organic_DD";

	// 1. Mesh
	vector< vector<double> > interface_nodes;
	vector<double> node(2);
	node[0] = 0.0;	node[1] = 1.0;	interface_nodes.push_back(node);
	node[1] = 0.0;	interface_nodes.push_back(node);

	double max_area = 0.01;
	my_mesh::MeshData mesh = my_mesh::generateMesh (filename, interface_nodes, max_area);


	// 2. Solve for solution
	arma::vec psi, n, p, u;
	double applied_psi = 0.0;
	solve_GummelIteration(mesh, psi, n, p, u, applied_psi);


	// 3. Plot solution
	bool is_to_plot=true;
	if (is_to_plot)
	{
		my_gnuplot::plot_ArmaVec(mesh, psi, filename+"_psi");
		my_gnuplot::plot_ArmaVec(mesh, n, filename+"_n");
		my_gnuplot::plot_ArmaVec(mesh, p, filename+"_p");
		my_gnuplot::plot_ArmaVec(mesh, u, filename+"_u");
	}

	return 0;
}
