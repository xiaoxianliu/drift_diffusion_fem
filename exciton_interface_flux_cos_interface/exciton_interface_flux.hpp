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

/* Functionals */
double interfaceIntegral(const my_mesh::MeshData& mesh, const arma::vec& u );

/* Interface flux */
int computeInterfaceNormalDerivative(	const my_mesh::MeshData &mesh, 
					const arma::vec &u,
					arma::vec &du_dnu_1, arma::vec &du_dnu_2);
//Note: Above, "nu" indicates the outer normal direction on surface



/* Plotting */
int plot_ArmaVec(const my_mesh::MeshData& mesh, const arma::vec u, const std::string& filename);
int plot_ArmaVec_on_Interface(const my_mesh::MeshData &mesh, const arma::vec &u, std::string filename);		

int plot_InterfaceArmaVec(const my_mesh::MeshData &mesh, const arma::vec &u, std::string filename);// Arma::vec defined on interface
int plot_InterfaceSTLVector(const my_mesh::MeshData &mesh, const std::vector<double> &u, std::string filename);
