#include <armadillo>
#include "mesh.hpp"


#ifndef linear_fem
#define linear_fem

namespace linear_fem
{
arma::mat assembleMatrixA(const my_mesh::MeshData& mesh, const arma::vec& a);
arma::mat assembleMatrixC(const my_mesh::MeshData& mesh, const arma::vec& c);
arma::mat assembleMatrixD(const my_mesh::MeshData& mesh, const arma::vec& d);
}


#endif
