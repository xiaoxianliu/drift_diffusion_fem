#include <iostream>
#include <armadillo>

#include "../triangle/mesh.hpp"
#include "fem_assemble.hpp"


//using namespace my_mesh;
//using namespace arma;

namespace linear_fem
{

arma::mat assembleMatrixGummel(const my_mesh::MeshData &mesh, const arma::vec &psi)
{
	int num_nodes = mesh.num_nodes;
	arma::mat M(num_nodes, num_nodes);
	M.zeros();




	return M;
}




}
