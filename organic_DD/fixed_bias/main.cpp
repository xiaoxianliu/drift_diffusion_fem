#include <iostream>
#include <vector>
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
/*	std::cout << "\n\n";
	std::cout << "lambda_squared is " << parameters::lambda_squared << "\n";
	std::cout << "epsilon_r = " << parameters::epsilon_rel << "\n";
	std::cout << "mu_n(donor) = " << parameters::mu_n_donor<<";\tmu_n(acceptor) = " << parameters::mu_n_acceptor 
		<< ";\tfield dependence parameter = " << parameters::gamma_n << "\n";
	std::cout << "mu_p(donor) = " << parameters::mu_p_donor<<"\tmu_p(acceptor) = " << parameters::mu_p_acceptor 
		<< ";\tfield dependence parameter = " << parameters::gamma_p << "\n";
	std::cout << "mu_x = " << parameters::mu_x << "\n";


*/



	string filename = "organic_DD";

	// 1. Mesh
	double x_offset = 0.0;
	double x_amplitude = -0.05;	// negative value means the interface varies in the "-x" direction
	int num_bumps = 0;
	double max_area = 2e-3;
	my_mesh::MeshData mesh = my_mesh::generateMesh_great_wall (	filename,
									x_offset,
									x_amplitude,
									num_bumps,
									max_area);






	// 2. Solve for solution
	// 2.1 solve
	arma::vec psi, n, p, u;
	double applied_psi = 0.0;
	solve_GummelIteration(mesh, psi, n, p, u, applied_psi);

	// 2.2 Plot solution
	bool is_to_plot=true;
	if (is_to_plot)
	{
		my_gnuplot::plot_ArmaVec(mesh, psi, filename+"_psi", "wxt");
		my_gnuplot::plot_ArmaVec(mesh, n, filename+"_n", "wxt");
		my_gnuplot::plot_ArmaVec(mesh, p, filename+"_p", "wxt");
		my_gnuplot::plot_ArmaVec(mesh, u, filename+"_u", "wxt");
//		my_gnuplot::plot_ArmaVec(mesh, p-n, filename+"_p-n", "wxt");

		my_gnuplot::plot_ArmaVec_on_DirichletBoundary(mesh, 1, psi, filename+"_psi(on_anode)");
		my_gnuplot::plot_ArmaVec_on_DirichletBoundary(mesh, 1, n, filename+"_n(on_anode)");
		my_gnuplot::plot_ArmaVec_on_DirichletBoundary(mesh, 1, p, filename+"_p(on_anode)");
		my_gnuplot::plot_ArmaVec_on_DirichletBoundary(mesh, 1, u, filename+"_u(on_anode)");
	}







	// 3. test
/*	bool is_to_test=false;
	if (is_to_test)
	{	arma::vec psi_test;
		test_Compute_psi(mesh, n, p, psi_test, applied_psi);
		std::cout << "difference between \"psi\" and \"psi_test\" is " << arma::norm(psi - psi_test, "inf") << "\n";
		my_gnuplot::plot_ArmaVec(mesh, psi_test, filename+"_psi_test", "wxt");
	}
*/




	// 4. Compute derived quantities and plot
	// 4.1 Electric field
	bool is_to_plot_E = false;
	if (is_to_plot_E)
	{	arma::vec E_x, E_y;
		compute_ElectricField(mesh, psi, E_x, E_y);
		my_gnuplot::plot_ArmaVec(mesh, E_x, filename+"_Ex", "png");
		my_gnuplot::plot_ArmaVec(mesh, E_y, filename+"_Ey", "png");
	}
	// 4.2 Fluxes
	bool is_to_plot_domain_flux = false;
	if (is_to_plot_domain_flux)
	{	// Flux of electrons
		arma::vec Fn_x, Fn_y;
		compute_Flux_n(mesh, n, psi, Fn_x, Fn_y);
		my_gnuplot::plot_ArmaVec(mesh, Fn_x, filename+"_Fn_x", "png");
		my_gnuplot::plot_ArmaVec(mesh, Fn_y, filename+"_Fn_y", "png");

		// Flux of holes
		arma::vec Fp_x, Fp_y;
		compute_Flux_p (mesh, p, psi, Fp_x, Fp_y);
		my_gnuplot::plot_ArmaVec(mesh, Fp_x, filename+"_Fp_x", "png");
		my_gnuplot::plot_ArmaVec(mesh, Fp_y, filename+"_Fp_y", "png");

		// Current = Fp - Fn
		my_gnuplot::plot_ArmaVec(mesh, Fp_x - Fn_x, filename+"_F_total_x", "png");

		// Flux of exciton
		arma::vec Fx_x, Fx_y;
		compute_Flux_x(mesh, u, Fx_x, Fx_y);
		my_gnuplot::plot_ArmaVec(mesh, Fx_x, filename+"_F_exciton", "png");
	}

	// 4.3 Boundary flux
	bool is_to_plot_boundary_flux = true;
	if (is_to_plot_boundary_flux)
	{	arma::vec Fn_nu1_anode, Fp_nu1_anode;
		compute_Boundary1Flux_n(mesh, n, psi, Fn_nu1_anode);
		compute_Boundary1Flux_p(mesh, p, psi, Fp_nu1_anode);
		my_gnuplot::plot_BoundaryArmaVec_on_DirichletBoundary(mesh, 1, Fn_nu1_anode, filename+"_Fn_nu1_anode", "png");
		my_gnuplot::plot_BoundaryArmaVec_on_DirichletBoundary(mesh, 1, Fp_nu1_anode, filename+"_Fp_nu1_anode", "png");
	}

	// 4.4 Boundary current density
	bool is_to_compute_boundary_current_density = true;
	if (is_to_compute_boundary_current_density)
	{	double Jp_anode, Jn_anode;
		Jn_anode = compute_Boundary1CurrentDensity_n(mesh, n, psi);
		Jp_anode = compute_Boundary1CurrentDensity_p(mesh, p, psi);

		std::cout << "Total current density on anode is J = " << Jp_anode + Jn_anode << "\n";
	}


return 0;

}
