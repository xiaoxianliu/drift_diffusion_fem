#include <string>
#include <armadillo>
#include "../triangle/mesh.hpp"


#ifndef LINEAR_FEM_H
#define LINEAR_FEM_H

namespace linear_fem
{

// Assemble FEM coefficient matrix for different differential operators (1st order approximation)
arma::mat assembleMatrixA(const my_mesh::MeshData& mesh, const arma::vec& a);
arma::mat assembleMatrixC(const my_mesh::MeshData& mesh, const arma::vec& c);
arma::mat assembleMatrixD(const my_mesh::MeshData& mesh, const arma::vec& d);
arma::mat assembleMatrixGummel(const my_mesh::MeshData& mesh, const arma::vec& psi, const arma::vec &mu);

// interpolation at nodal points of a given mesh
arma::vec interpolateConstant(const my_mesh::MeshData &mesh, double C);		// interpolation of a constant
arma::vec interpolateFunction(const my_mesh::MeshData& mesh, double (*f)(double, double));	// interpolation of a function

// L2 projection of a given function/vector
arma::vec L2project_Func(	const my_mesh::MeshData &mesh, double (*f)(double, double) );
arma::vec L2project_Vec(	const my_mesh::MeshData &mesh, const arma::vec &u);


// Functionals
double integrate_Domain(const my_mesh::MeshData &mesh, const arma::vec u);

}


#endif
