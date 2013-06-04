#include <string>

#include <armadillo>
#include "../triangle/mesh.hpp"


#ifndef EXCITON_H
#define EXCITON_H

/* Mesh */
my_mesh::MeshData generateMesh(std::string filename);

/* Functions associated with exciton equation */
double generation_rate(double x, double y);			// exciton generation
double diffusion_coefficient(double x, double y);		// exciton diffusion coefficient
double decay_rate(double x, double y);				// exciton decay rate (= 1/ life_time)
double interface_reaction_rate(double x, double y);		// exciton reaction rate on interface

double uD_left(double x, double y);			// Dirichlet boundary condition on the left boundary
double uD_right(double x, double y);			// Dirichlet boundary condition on the right boundary


/* Linear interpolation of functions */
arma::vec interpolate_func(const my_mesh::MeshData &mesh, double (*f)(double, double));

/* Linear projection of functions */
arma::vec project_func(const my_mesh::MeshData &mesh, double (*f)(double, double));



#endif
