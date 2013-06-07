#include <string>
#include <vector>

#include <armadillo>

#include "../triangle/mesh.hpp"

/* Mesh */
my_mesh::MeshData generateMesh_CosInterface(const std::string &filename, double y_control);


/* Function approximation */
/* Function interpolation at mesh nodes */
arma::vec interpolateFunc(const my_mesh::MeshData& mesh, double (*f)(double, double));


/* State equation */
arma::vec solveStateEq(const my_mesh::MeshData mesh);
/* Adjoint equation */
arma::vec solveAdjointEq(const my_mesh::MeshData &mesh);



/* Shape gradient */
arma::vec computeShapeGradient(	const my_mesh::MeshData &mesh, const arma::vec u, const arma::vec xi);



/* Plotting */
int plot_ArmaVec(const my_mesh::MeshData& mesh, const arma::vec u, const std::string& filename);
int plot_ArmaVec_Interface(const my_mesh::MeshData &mesh, const arma::vec &u, std::string filename);		
int plot_STLVector_Interface(const my_mesh::MeshData &mesh, const std::vector<double> &u, std::string filename);
