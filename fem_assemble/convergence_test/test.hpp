#include <string>
#include "../../triangle/mesh.hpp"
#include <armadillo>

//my_mesh::MeshData generateMesh_CosInterface(const std::string &filename, double y_control, double max_area);
int solveEq(const my_mesh::MeshData mesh, const std::string &filename, arma::vec &u, double &error);
