#include <string>
#include <armadillo>

#include "../triangle/mesh.hpp"

namespace my_gnuplot
{

// Plot arma::vec over the entire domain "Omega" 
int plot_ArmaVec(const my_mesh::MeshData& mesh, const arma::vec u, const std::string& filename, std::string terminal="png");

/**********************************************************************************************************************************/
// Plot on interface
// (1) Plot a arma::vec defined over "Omega" only on the interface "Gamma" 
int plot_ArmaVec_on_Interface(const my_mesh::MeshData &mesh, const arma::vec &u, std::string filename, std::string terminal="png");


// (2) Plot arma::vec which is only defined on the interface Gamma 
int plot_InterfaceArmaVec(const my_mesh::MeshData &mesh, const arma::vec &u, std::string filename, std::string terminal="png");

// (3) Plot std::vector which is only defined on the interface Gamma 
int plot_InterfaceSTLVector(const my_mesh::MeshData &mesh, const std::vector<double> &u, std::string filename, std::string terminal="png");

/**********************************************************************************************************************************/
// Plot on Dirichlet boundary
int plot_ArmaVec_on_DirichletBoundary(	const my_mesh::MeshData &mesh,
					int boundary_marker,
					const arma::vec &u, 
					std::string filename, 
					std::string terminal = "png");

int plot_BoundaryArmaVec_on_DirichletBoundary(	const my_mesh::MeshData &mesh,
						int boundary_marker, 
						const arma::vec &u, 
						std::string filename, 
						std::string terminal="png");

}
