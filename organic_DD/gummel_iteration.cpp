#include <iostream>
#include <cmath>

#include <armadillo>

#include "../triangle/mesh.hpp"
#include "parameters.hpp"
#include "gummel_iteration.hpp"





/*************************************************************************************************************************/
// Declaration of local functions
// initialize vectors with a guess of solution to Drift-Diffusion equations for organic semiconductor
int initialization(	const my_mesh::MeshData &mesh,
			arma::vec &psi,
			arma::vec &n,
			arma::vec &p,
			arma::vec &u,
			double applied_psi);



/*************************************************************************************************************************/
// Main function
int solve_GummelIteration(	const my_mesh::MeshData &mesh,
				arma::vec &psi,
				arma::vec &n,
				arma::vec &p,
				arma::vec &u,
				double applied_psi			// applied potential at "anode"
			)
{
	// 1. Initialize {psi, n, p, u} on given mesh
	initialization(mesh, psi, n, p, u, applied_psi);

	// 2. Gummel's iteration
	// 2.1 setup
	arma::vec prev_psi = psi, prev_n = n, prev_p = p, prev_u = u;
	arma::vec new_psi, new_n, new_p, new_u;

	double gummel_err = 1.0, gummel_tol = 1e-5;
	int max_iter = 0;

	// 2.2 gummel's iteration
	for (int iter=0; iter<max_iter; iter++)
	{
		std::cout << "Gummel iteration: " << iter << "\n";

		// 2.2.1 Solve equations one by one
		solve_ContinuityEq_n(mesh, prev_psi, prev_n, prev_p, prev_u, new_n);
		solve_ContinuityEq_p(mesh, prev_psi, prev_n, prev_p, prev_u, new_p);
		solve_ContinuityEq_x(mesh, prev_psi, prev_n, prev_p, prev_u, new_u);
		solve_NonlinearPoissonEq(mesh, prev_psi, n, p, u, new_psi, applied_psi);

		// 2.2.2 compute gummel error
		double max_dn = arma::norm(new_n-prev_n, "inf");
		double max_dp = arma::norm(new_p-prev_p, "inf");
		double max_du = arma::norm(new_u-prev_u, "inf");
		double max_dpsi = arma::norm(new_psi-prev_psi, "inf");

		double max_dn_vs_dp = fmax(max_dn, max_dp);
		double max_du_vs_dpsi = fmax(max_du, max_dpsi);
		gummel_err = fmax(max_dn_vs_dp, max_du_vs_dpsi);

		std::cout << "Gummel's error is " << gummel_err << "\n";

		if (gummel_err < gummel_tol)
			break;

		prev_psi = new_psi;	prev_n = new_n;	prev_p = new_p;	prev_u = new_u;
	}

	psi = new_psi;	n = new_n;	p = new_p;	u = new_u;


	return 0;
}







/*************************************************************************************************************************/

// Definition of local functions

// Initialize solution vectors with a guess to start Gummel's iteration
int initialization(	const my_mesh::MeshData &mesh,
			arma::vec &psi,
			arma::vec &n,
			arma::vec &p,
			arma::vec &u,
			double applied_psi)
{
/**** Initialize every unkown with linear interpolation of the boundary values ************/
using namespace parameters;

	int num_nodes = mesh.num_nodes;
	psi.resize(num_nodes);	n.resize(num_nodes);	p.resize(num_nodes);	u.resize(num_nodes);
	psi.zeros();	n = n.zeros();	p.zeros();	u.zeros();

	for (int i=0; i<mesh.num_nodes; i++)
	{	double x = mesh.nodes[i][0];

		psi(i) = psi_A + applied_psi + (psi_C - psi_A - applied_psi)*(x+1)/2.0;
		n(i) = n_A + (n_C - n_A)*(x+1)/2.0;
		p(i) = p_A + (p_C - p_A)*(x+1)/2.0;
		u(i) = u_A + (u_C - u_A)*(x+1)/2.0;
	}

	return 0;
}
