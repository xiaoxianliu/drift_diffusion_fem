#include <string>
#include <armadillo>
#include "../triangle/mesh.hpp"


#ifndef LINEAR_FEM_H
#define LINEAR_FEM_H

namespace linear_fem
{
arma::mat assembleMatrixA(const my_mesh::MeshData& mesh, const arma::vec& a);
arma::mat assembleMatrixC(const my_mesh::MeshData& mesh, const arma::vec& c);
arma::mat assembleMatrixD(const my_mesh::MeshData& mesh, const arma::vec& d);

// interpolate function at nodal points of given mesh
arma::vec interpolateFunction(const my_mesh::MeshData& mesh, double (*f)(double, double));	

}


#endif
