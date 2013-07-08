#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>

#include <armadillo>

#include "../../triangle/mesh.hpp"
#include "../../gnuplot/my_gnuplot.hpp"


#include "gummel_iteration.hpp"
#include "auxillary.hpp"
#include "test.hpp"
#include "flux_and_gradient.hpp"


int main (int argc, char* argv[])
{
using namespace std;

	// 0. Output some parameters
	std::cout << "\n\n";
	std::cout << "lambda_squared is " << parameters::lambda_squared << "\n";
	std::cout << "epsilon_r = " << parameters::epsilon_rel << "\n";
	std::cout << "mu_n(donor) = " << parameters::mu_n_donor<<";\tmu_n(acceptor) = " << parameters::mu_n_acceptor 
		<< ";\tfield dependence parameter = " << parameters::gamma_n << "\n";
	std::cout << "mu_p(donor) = " << parameters::mu_p_donor<<"\tmu_p(acceptor) = " << parameters::mu_p_acceptor 
		<< ";\tfield dependence parameter = " << parameters::gamma_p << "\n";
	std::cout << "mu_x = " << parameters::mu_x << "\n";

	string filename = "organic_DD";






	// 1. Mesh
	vector< vector<double> > interface_nodes;
	vector<double> node(2);
	node[0] = 0.0;	node[1] = 1.0;	interface_nodes.push_back(node);
	node[1] = 0.0;	interface_nodes.push_back(node);

	double max_area = 2e-3;
	my_mesh::MeshData mesh = my_mesh::generateMesh (filename, interface_nodes, max_area);




	// 2. Compute I-V characteristics
	// 2.1 Various applied potential
	int num_applied_psi = 20;
	double delta_psi = 2.0;
	arma::vec applied_psi_s (num_applied_psi);
	for (int i=0; i<applied_psi_s.n_rows; i++){
		applied_psi_s(i) = i * delta_psi;
	}
	arma::vec Js(num_applied_psi);	Js.zeros();

	// 2.2 solve for current for each applied_psi
	for (int i=0; i<applied_psi_s.n_rows; i++)
	{
		double applied_psi = applied_psi_s(i);
		std::cout << "Applied psi #" << i << ": " << applied_psi << "\n";

		// 2.2.1 Solve for densities and potentials
		arma::vec psi, n, p, u;
		solve_GummelIteration(mesh, psi, n, p, u, applied_psi);

		// 2.2.2 Compute current density on either "anode" or "cathode"

		double Jp_anode, Jn_anode;
		Jn_anode = compute_Boundary1CurrentDensity_n(mesh, n, psi);
		Jp_anode = compute_Boundary1CurrentDensity_p(mesh, p, psi);

		Js(i) = Jn_anode + Jp_anode;		
	}

	std::cout << "Js = \n" << Js << "\n";

	// 2.3 Plot I-V characteristics
	// 2.3.1 Data
	string iv_dat_filename = filename + "_iv.dat";
	ofstream iv_dat_stream;
	iv_dat_stream.open (iv_dat_filename.c_str());
	for (int i=0; i<Js.n_rows; i++){
		iv_dat_stream << applied_psi_s(i) << "\t" << Js(i) << "\n";
	}
	iv_dat_stream.close();

	// 2.3.2 write gnuplot file
	string iv_gnuplot_filename = filename + "_iv.gnuplot";
	ofstream iv_gnuplot_stream;
	iv_gnuplot_stream.open(iv_gnuplot_filename.c_str());
	iv_gnuplot_stream << "set terminal png\n";
	iv_gnuplot_stream << "set output \"" << filename + "_iv.png\"\n";
	iv_gnuplot_stream << "plot \"" << iv_dat_filename << "\" with lines\n";
	iv_gnuplot_stream.close();

	// 2.3.3 run "gnuplot" from command line
	string cmd = "gnuplot " + iv_gnuplot_filename + " --persist\n";
	system(cmd.c_str());
	




return 0;

}
