#include <vector>

#include <armadillo>
#include "../triangle/mesh.hpp"

double diffusion_coefficient(double x, double y)
{	return 1.0;
}

double decay_rate(double x, double y)
{	return 1.0;
}

double generation_rate(double x, double y)
{	return 1.0;
}

double interface_reaction_rate(double x, double y)
{	return 1.0e3;
}


/* boundary conditions */
double uD_left(double x, double y)
{
	return 0.0;
}

double uD_right(double x, double y)
{	return 0.0;
}




/* Interpolate functions */
arma::vec interpolate_func(const my_mesh::MeshData &mesh, double (*f)(double, double))
{
	arma::vec u(mesh.num_nodes);
	u.zeros();

	for (int v=0; v<mesh.num_nodes; v++)
	{	double x = mesh.nodes[v][0];
		double y = mesh.nodes[v][1];

		u(v) = f(x,y);
	}

	return u;
}


/* Projection of functions to linear FEM space */
/* int_Omega ( phi_test, f) *dx, where "phi_test" is a basis function in FEM space */
arma::vec project_func(const my_mesh::MeshData &mesh, double (*f)(double, double))
{
	arma::vec u(mesh.num_nodes);
	u.zeros();

	for (int v=0; v<mesh.num_nodes; v++)
	{	double x = mesh.nodes[v][0];
		double y = mesh.nodes[v][1];

		std::vector<int> Ts = mesh.topology0to2[v];
		for (int i=0; i<Ts.size(); i++)
		{	int t = Ts[i];				// element index "t"
			double area_T = mesh.ele_areas[t];

			u(v) += f(x,y) * area_T / 3.0;
		}

	}
	return u;
}






