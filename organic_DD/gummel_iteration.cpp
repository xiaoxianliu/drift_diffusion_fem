#include <iostream>
#include <armadillo>

#include "../triangle/mesh.hpp"
#include "parameters.hpp"
#include "main.hpp"


int initialization(	const my_mesh::MeshData &mesh,
			arma::vec &psi,
			arma::vec &n,
			arma::vec &p,
			arma::vec &u,
			double applied_psi);





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
	arma::vec prev_psi = psi,
		  prev_n = n,
		  prev_p = p,
		  prev_u = u;
	double gummel_err = 1.0, gummel_tol = 1e-3;
	int max_iter = 5;
	for (int iter=0; iter<max_iter; iter++)
	{

//		given "psi", 
//(1)compute "-E = grad(psi)"
//(2)compute normal derivative of "psi" on interface


//		solve_NContinuityEq();
//		solve_PContinuityEq();
//		solve_XContinuityEq();
//		solve_NonlinearPoissonEq();

//		compute gummel error

		if (gummel_err < gummel_tol)
			break;
	}


	return 0;
}









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
