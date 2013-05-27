#include <armadillo>
#include "mesh.hpp"

using namespace arma;

/* Mesh */
int polyfileSmartInterface(const string polyname, const int num_squares);
int TriangleMesh(MeshData &mesh, string& filename);

/* Linear system */
// assembling linear system
mat assemble_matrix_A(const MeshData &mesh, const vec &a);			// A_ij = a*( grad(phi_i), grad(phi_j) )*dx(Omega)
mat assemble_matrix_C(const MeshData& mesh, const arma::vec &c);		// C_ij = c* phi_i * phi_j *dx(Omega)
mat assemble_matrix_D(const MeshData& mesh, const arma::vec &d);		// D_ij = d* phi_i * phi_j *ds(interface)
vec assemble_rhsvector(const MeshData& mesh, const vec &func);

// Boundary condition
void applyDirichletBC(const MeshData &mesh, arma::mat &A, arma::vec &rhs);

// Interpolation functions
vec interpolateCoefficient_a(const MeshData &mesh);
vec interpolateCoefficient_c(const MeshData &mesh);
vec interpolateCoefficient_d(const MeshData &mesh);
vec interpolateRHS(const MeshData &mesh);

// plot solution
int PlotSolution(const MeshData& mesh, const arma::vec u, const string& polyname);

// Solve PDE by FEM
int Solve(string polyname, const int num_squares, double& interface_drainage);

// compute various functionals
double interfaceIntegral(const MeshData& mesh, const arma::vec& u, const arma::vec& d);

