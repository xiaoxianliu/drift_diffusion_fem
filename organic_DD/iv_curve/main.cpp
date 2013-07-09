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

	string filename = "organic_DD";


	// 1. Meshes with cosine interfaces of different period numbers
	double x_offset = 0.0;
	double x_amplitude = -0.4;	// negative value means the interface varies in the "-x" direction
	int num_interface_periods;
	int init_num_interface_nodes = 401;
	double max_area = 5e-3;

	my_mesh::MeshData mesh0, mesh1, mesh2;

	num_interface_periods = 0;	string filename0 = filename + "_0_period";
	mesh0 = my_mesh::generateMesh_cosine_interface (filename+"_0_period", x_offset, x_amplitude, num_interface_periods, 
							init_num_interface_nodes, max_area);
	num_interface_periods = 1;	string filename1 = filename + "_1_period";
	mesh1 = my_mesh::generateMesh_cosine_interface (filename1, x_offset, x_amplitude, num_interface_periods, 
							init_num_interface_nodes, max_area);
	num_interface_periods = 2;	string filename2 = filename + "_2_period";
	mesh2 = my_mesh::generateMesh_cosine_interface (filename2, x_offset, x_amplitude, num_interface_periods, 
							init_num_interface_nodes, max_area);




	// 2. Compute I-V characteristics
	// 2.1 Various applied potential
	int num_applied_psi = 20;
	double delta_psi = 2.0;
	arma::vec applied_psi_s (num_applied_psi);
	for (int i=0; i<applied_psi_s.n_rows; i++){
		applied_psi_s(i) = i * delta_psi;
	}

	arma::vec Js0(num_applied_psi);	Js0.zeros();
	arma::vec Js1(num_applied_psi);	Js1.zeros();
	arma::vec Js2(num_applied_psi);	Js2.zeros();

	// 2.2 solve for current for each applied_psi
	for (int i=0; i<applied_psi_s.n_rows; i++)
	{
		double applied_psi = applied_psi_s(i);
		std::cout << "Applied psi #" << i << ": " << applied_psi << "\n";

		arma::vec psi, n, p, u;
		double Jp_anode, Jn_anode;

		// 2.2.0 Mesh 0
		std::cout << "Computing mesh with flat interface\n";
		// 2.2.0.1 Solve for densities and potentials
		solve_GummelIteration(mesh0, psi, n, p, u, applied_psi);

		// 2.2.0.2 Compute current density on either "anode" or "cathode"
		Jn_anode = compute_Boundary1CurrentDensity_n(mesh0, n, psi);
		Jp_anode = compute_Boundary1CurrentDensity_p(mesh0, p, psi);
		Js0(i) = Jn_anode + Jp_anode;


		// 2.2.1 "Mesh 1"
		std::cout << "Computing mesh with cosine interface of 1 period\n";
		// 2.2.1.1 Solve for densities and potentials
		solve_GummelIteration(mesh1, psi, n, p, u, applied_psi);

		// 2.2.0.2 Compute current density on either "anode" or "cathode"
		Jn_anode = compute_Boundary1CurrentDensity_n(mesh1, n, psi);
		Jp_anode = compute_Boundary1CurrentDensity_p(mesh1, p, psi);
		Js1(i) = Jn_anode + Jp_anode;


		// 2.2.2 "Mesh 2"
		std::cout << "Computing mesh with cosine interface of 2 period\n";
		// 2.2.2.1 Solve for densities and potentials
		solve_GummelIteration(mesh2, psi, n, p, u, applied_psi);

		// 2.2.2.2 Compute current density on either "anode" or "cathode"
		Jn_anode = compute_Boundary1CurrentDensity_n(mesh2, n, psi);
		Jp_anode = compute_Boundary1CurrentDensity_p(mesh2, p, psi);
		Js2(i) = Jn_anode + Jp_anode;

		std::cout << "\n\n";
	}


	// 2.3 Plot I-V characteristics
	// 2.3.1 Data
	string iv_dat_filename = filename + "_iv.dat";
	ofstream iv_dat_stream;
	iv_dat_stream.open (iv_dat_filename.c_str());
	for (int i=0; i<num_applied_psi; i++){
		iv_dat_stream 
			<< applied_psi_s(i) << "\t" 
			<< Js0(i) << "\t" 
			<< Js1(i) << "\t" 
			<< Js2(i) <<"\n";
	}
	iv_dat_stream.close();

	// 2.3.2 write gnuplot file
	string iv_gnuplot_filename = filename + "_iv.gnuplot";
	ofstream iv_gnuplot_stream;
	iv_gnuplot_stream.open(iv_gnuplot_filename.c_str());
	iv_gnuplot_stream << "set terminal png\n";
	iv_gnuplot_stream << "set output \"" << filename + "_iv.png\"\n";
	iv_gnuplot_stream << "plot\t"
				<< "\"" << iv_dat_filename << "\" u 1:2 t \"" << filename0 << "\"" << " w linespoints,\\\n" 
				<< "\"" << iv_dat_filename << "\" u 1:3 t \"" << filename1 << "\"" << " w linespoints,\\\n" 
				<< "\"" << iv_dat_filename << "\" u 1:4 t \"" << filename2 << "\"" << " w linespoints" 
				<< "\n";
	iv_gnuplot_stream.close();

	// 2.3.3 run "gnuplot" from command line
	string cmd = "gnuplot " + iv_gnuplot_filename + " --persist\n";
	system(cmd.c_str());
	




return 0;

}
