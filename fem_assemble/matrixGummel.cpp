#include <iostream>
#include <cmath>
#include <armadillo>

#include "../triangle/mesh.hpp"
#include "fem_assemble.hpp"


//using namespace my_mesh;
//using namespace arma;

double BernoulliFunc(double x);

namespace linear_fem
{

arma::mat assembleMatrixGummel(const my_mesh::MeshData &mesh, const arma::vec &psi, const arma::vec &mu )
{
	int num_nodes = mesh.num_nodes;
	arma::mat M(num_nodes, num_nodes);
	M.zeros();

	for (int t=0; t<mesh.elements.size(); t++)			// element "t"
	{
		double mu_t = mu(t);					// element-wise constant mobility "mu_t"
		double area_t = mesh.ele_areas[t];			// area of triangle "t"

		for (int i=0; i<3; i++)					// cycle through all 3 vertices as "row index"
		{							// compute M(v0, v0), M(v0, v1), M(v0, v2)
			int v0 = mesh.elements[t][i];
			int v1 = mesh.elements[t][(i+1)%3];
			int v2 = mesh.elements[t][(i+2)%3];

			double x0 = mesh.nodes[v0][0], y0 = mesh.nodes[v0][1],
				x1 = mesh.nodes[v1][0], y1 = mesh.nodes[v1][1],
				x2 = mesh.nodes[v2][0], y2 = mesh.nodes[v2][1];

			double psi0 = psi(v0), psi1 = psi(v1), psi2 = psi(v2);


			arma::mat J(2,2);				// Jacobian matrix "d_x/d_xi" of transformation
			J(0,0) = x1 - x0;	J(1,0) = y1 - y0;
			J(0,1) = x2 - x0;	J(1,1) = y2 - y0;

			arma::mat JtJ_inv(2,2);				// (J_transpose * J) ^(-1)
			JtJ_inv = ( J.t() * J ).i();

			arma::vec dphi_dxi = - arma::ones<arma::vec>(2);	// "dphi_dxi = (-1, -1)^t"
			arma::vec b0(2), b1(2), b2(2);
			b0(0) = -BernoulliFunc(psi1 - psi0);	b0(1) = -BernoulliFunc(psi2 - psi0);
			b1(0) = BernoulliFunc(psi1 - psi0) * std::exp(psi1 - psi0);	b1(1) = 0.;
			b2(0) = 0.;	b2(1) = BernoulliFunc(psi2 - psi0) * std::exp(psi2 - psi0);

			M(v0, v0) += mu_t * area_t * arma::dot( b0, JtJ_inv * dphi_dxi );
			M(v0, v1) += mu_t * area_t * arma::dot( b1, JtJ_inv * dphi_dxi );
			M(v0, v2) += mu_t * area_t * arma::dot( b2, JtJ_inv * dphi_dxi );
		}
		

	}


	return M;
}




}




double BernoulliFunc(double x)
{
	if ( x < 1e-10 && x>-1e-10 )
		return 1.0;
	else
	{
		double denominator = std::exp(x) - 1.0;
		return x/denominator;
	}
}
