#include <string>
#include <vector>

#include <armadillo>

#include "../triangle/mesh.hpp"

/* Mesh */
my_mesh::MeshData generateMesh_CosInterface(const std::string &filename, double y_control);
/* interface curvature */
std::vector<double> interfaceCurvature(const my_mesh::MeshData &mesh);


/* Function approximation */
/* Function interpolation at mesh nodes */
arma::vec interpolateFunc(const my_mesh::MeshData& mesh, double (*f)(double, double));


/* State equation */
arma::vec solveStateEq(const my_mesh::MeshData mesh);
/* Adjoint equation */
arma::vec solveAdjointEq(const my_mesh::MeshData &mesh);



/* Plotting */
int plotSolutionVec(const my_mesh::MeshData& mesh, const arma::vec u, const std::string& filename);

