#include <iostream>
#include <cstdlib>
#include <fstream>
#include <armadillo>


#include "mesh.hpp"
#include "exciton_diffusion.hpp"

//using namespace std;
using namespace arma;

/********** Declaration of auxillary functions *****************************/
/* defined at the bottom of this file */
double aFunc(double x, double y);				// a(x) in "-div( a*grad(u) ) + c*u = f"
double cFunc(double x, double y);				// c(x) in "-div( a*grad(u) ) + c*u = f"
double InterfaceRateFunc(double x, double y);			// reaction rate on interface
double RHSFunc(double x, double y);				// f(x) in "-div( a*grad(u) ) + c*u = f"
double DirichletFunc(double x, double y);			// u_D(x) on Dirichlet boundaries




/********** Definitions of functions starts ********************************/

/* 1. Assemble coefficient matrix */

/* 1.1 "Laplacian" operator:
	i.e. A_ij = (grad_phi_i, grad_phi_j) *dx
Note: no boundary conditions are imposed yet.
*/
mat assemble_matrix_A(const MeshData &mesh, const vec &a)
{
	int num_nodes = mesh.num_nodes;

	/* Initialize matrix to all 0's*/
	mat A=arma::zeros(num_nodes, num_nodes);

	/* Define matrix entries */
	for (int t=0; t<mesh.num_elements; t++)			// loop through all elments
	{	if (mesh.num_nodes_per_ele!=3)			// 3 for linear element; 6 for quadratic element, etc
			{std::cout << "Element has to be linear... i.e. 3 nodes per element " << std::endl; exit(1);}
		/* Vertices of triangle "t" */
		int v0 = mesh.elements[3*t];
		int v1 = mesh.elements[3*t+1];
		int v2 = mesh.elements[3*t+2];

		vec r0(2), r1(2), r2(2);			// 2-dim vectors of coordinates;
		r0(0) = mesh.nodes[2*v0];	r0(1) = mesh.nodes[2*v0+1];
		r1(0) = mesh.nodes[2*v1];	r1(1) = mesh.nodes[2*v1+1];
		r2(0) = mesh.nodes[2*v2];	r2(1) = mesh.nodes[2*v2+1];

		/* factor shared by all "A_ij(t)" */
		double factor;
		factor = (a(v0) + a(v1) + a(v2))/3.0 /(4.0*mesh.ele_areas[t]);

		/* Modify the 9 entries corresponding to the 3 vertices */
		A(v0,v0) += factor * dot(r1-r2, r1-r2);
		A(v0,v1) += factor * (- dot(r0-r2, r1-r2));
		A(v0,v2) += factor * (- dot(r2-r1, r0-r1));

		A(v1,v1) += factor * dot(r0-r2, r0-r2);
		A(v1,v0) += factor * (- dot(r0-r2, r1-r2));
		A(v1,v2) += factor * (- dot(r1-r0, r2-r0));

		A(v2,v2) += factor * dot(r0-r1, r0-r1);
		A(v2,v1) += factor * (- dot(r1-r0, r2-r0));
		A(v2,v0) += factor * (- dot(r2-r1, r0-r1));

	}

	return A;
}

/* 1.2 C_ij = c* phi_i * phi_j *dx(Omega) */
mat assemble_matrix_C(const MeshData& mesh, const arma::vec &c)
{
	/* Initialize C to a zero matrix of the right size */
	arma::mat C = arma::zeros(mesh.num_nodes, mesh.num_nodes);

	/* Define matrix entries */
	for (int t=0; t<mesh.num_elements; t++)
	{
		int v0 = mesh.elements[3*t];
		int v1 = mesh.elements[3*t+1];
		int v2 = mesh.elements[3*t+2];

		C(v0, v0) += mesh.ele_areas[t] * c(v0) / 3.0;
		C(v1, v1) += mesh.ele_areas[t] * c(v1) / 3.0;
		C(v2, v2) += mesh.ele_areas[t] * c(v2) / 3.0;
	}

	return C;	 
}

/* 1.3 D_ij = d* phi_i * phi_j *ds(interface) (related to interface) */
mat assemble_matrix_D(const MeshData& mesh, const arma::vec &d)
{
	arma::mat D(mesh.num_nodes, mesh.num_nodes);
	D.zeros();
	for (int e=0; e<mesh.num_edges; e++)
	{	if (mesh.edge_markers[e]==5)
		{
			int v0 = mesh.edges[2*e];
			int v1 = mesh.edges[2*e+1];
			double length = mesh.edge_lengths[e];

			D(v0, v0) = d(v0) * length / 2.0;
			D(v1, v1) = d(v1) * length / 2.0;
		}
	}
	return D;
}


/* 2. Assemble right-hand side vector for linear system */
vec assemble_rhsvector(const MeshData& mesh, const vec &func)
{	vec rhs(mesh.num_nodes);
	rhs.zeros();

	for (int t=0; t<mesh.num_elements; t++)
	{	if (mesh.num_nodes_per_ele!=3)			// 3 for linear element; 6 for quadratic element, etc
			{std::cout << "Element has to be linear... i.e. 3 nodes per element " << std::endl; exit(1);}
	/* Vertices of triangle "t" */
		int v0 = mesh.elements[3*t];
		int v1 = mesh.elements[3*t+1];
		int v2 = mesh.elements[3*t+2];

		rhs(v0) += mesh.ele_areas[t] * func(v0) / 3.0;
		rhs(v1) += mesh.ele_areas[t] * func(v1) / 3.0;
		rhs(v2) += mesh.ele_areas[t] * func(v2) / 3.0;
	}

	/* Apply Dirichlet boundary condition  */

	return rhs;
}


/* 3. Dirichlet boundary conditions */
void applyDirichletBC(const MeshData &mesh, arma::mat &M, arma::vec &rhs)
{
	/* Modify coefficient matrix: change rows of Dirichlet nodes v_i to (0,0,..., 1, ...0,0) with "1" at v_i-th entry */
	for (int e=0; e<mesh.num_edges; e++)
	{	int marker = mesh.edge_markers[e];

		/* by default, 1 and 3 coorespond to "left" and "right" both of which are Dirichlet boundaries. */
		if (marker==1 || marker==3)			
		{	int v0 = mesh.edges[2*e];
			int v1 = mesh.edges[2*e+1];
			M(v0, span::all) = zeros<mat>(1, mesh.num_nodes);
			M(v0,v0) = 1.0;
			M(v1, span::all) = zeros<mat>(1, mesh.num_nodes); 
			M(v1,v1) = 1.0;
		}
	}

	/* Modify right-hand side vector: change entries of Dirichlet nodes v_i to its Dirichlet boundary value */
	for (int e=0; e<mesh.num_edges; e++)
	{	int marker = mesh.edge_markers[e];
		if (marker==1 || marker==3)
		{	double x, y;

			int v0 = mesh.edges[2*e];
			x = mesh.nodes[2*v0];
			y = mesh.nodes[2*v0+1];
			rhs(v0) = DirichletFunc(x,y);					// preassumed Dirichlet boundary values

			int v1 = mesh.edges[2*e+1];
			x = mesh.nodes[2*v1];
			y = mesh.nodes[2*v1+1];
			rhs(v1) = DirichletFunc(x,y);					// preassumed Dirichlet boundary values
		}
	}

}


/* 4. Interpolated coefficient and functions */
// 4.1 linearly interpolated a(x) in "-div( a*grad(u) ) + c*u = f"
vec interpolateCoefficient_a(const MeshData &mesh)
{	vec a_vec(mesh.num_nodes);
	for (int i=0; i<mesh.num_nodes; i++)
	{	double x = mesh.nodes[2*i];
		double y = mesh.nodes[2*i+1];
		a_vec(i) = aFunc(x,y);
	}
	return a_vec;	
}

// 4.2 linearly interpolated c(x) in "-div( a*grad(u) ) + c*u = f"
vec interpolateCoefficient_c(const MeshData &mesh)
{	vec c_vec(mesh.num_nodes);
	for (int i=0; i<mesh.num_nodes; i++)
	{	double x = mesh.nodes[2*i];
		double y = mesh.nodes[2*i+1];
		c_vec(i) = cFunc(x,y);
	}
	return c_vec;	
}

// 4.3 linearly interpolated "reaction rate" on left-right interface"
vec interpolateCoefficient_d(const MeshData &mesh)		// Contribution of interface reaction to the coefficient matrix
{	vec d_vec(mesh.num_nodes);
	d_vec.zeros();
	for (int e=0; e<mesh.num_edges; e++)
	{	if (mesh.edge_markers[e]==5)			// default value of edge marker for interface
		{	double x, y;

			int v0 = mesh.edges[2*e];
			x = mesh.nodes[2*v0];	y = mesh.nodes[2*v0+1];
			d_vec(v0) = InterfaceRateFunc(x,y);

			int v1 = mesh.edges[2*e+1];
			x = mesh.nodes[2*v1];	y = mesh.nodes[2*v1+1];
			d_vec(v1) = InterfaceRateFunc(x,y);
		}
	}
	return d_vec;	
}

// 4.4 linearly interpolated f(x) in "-div( a*grad(u) ) + c*u = f"
vec interpolateRHS(const MeshData &mesh)
{	vec f_vec(mesh.num_nodes);
	for (int i=0; i<mesh.num_nodes; i++)
	{	double x = mesh.nodes[2*i];
		double y = mesh.nodes[2*i+1];
		f_vec(i) = RHSFunc(x,y);
	}
	return f_vec;
}



/**************** Below are the functions only used within this file *************************************/
/* 5. Auxillary functions of original equations, not included in "exciton_diffusion.hpp" for declaration */
double aFunc(double x, double y)				// a(x) in "-div( a*grad(u) ) + c*u = f"
{	return 1.0;
}

double cFunc(double x, double y)				// c(x) in "-div( a*grad(u) ) + c*u = f"
{	return 1;
}

double InterfaceRateFunc(double x, double y)
{	return 1e1;
}
double RHSFunc(double x, double y)				// f(x) in "-div( a*grad(u) ) + c*u = f"
{
//	double f = -1.0 + 0.5*x*x;			// corresponds to u=0.5*x^2
//	double f = exp(-x*x) * (-4*x*x +3);		// corresponds to u=exp(-x^2)
//	double f = x>=0 ? (x+2):(-x+2);			//corresponds to u = |x|+2
	double f = 1.0;					// homogeneous exciton photo-generation (light)
	return f;
}
double DirichletFunc(double x, double y)			// u_D(x) on Dirichlet boundaries
{
//	double u_D = 0.5*x*x;				// cooresponds to u=0.5*x^2
//	double u_D = exp(-x*x);				// corresponds to u = exp(-x^2)
//	double u_D = x>=0? (x+2):(-x+2);		// corresponds to u = |x| + 2
	double u_D = 0.0;				// homogeneous boundary condition
	return u_D;
}

