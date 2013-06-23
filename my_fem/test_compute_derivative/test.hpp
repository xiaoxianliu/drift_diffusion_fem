#include <string>
#include <armadillo>
#include "../../triangle/mesh.hpp"
#include "../my_fem.hpp"



// Functions for test
double f(double x, double y);
double grad_f_x(double x, double y);
double grad_f_y(double x, double y);

// Function that computes gradients of given function
int computeGradient(	const my_mesh::MeshData &mesh,
			const arma::vec &psi,
			arma::vec &psi_grad_x,
			arma::vec &psi_grad_y);


// Plotting functions
/* Plot arma::vec over the entire domain "Omega" */
int plot_ArmaVec (const my_mesh::MeshData& mesh, const arma::vec u, const std::string& filename);
