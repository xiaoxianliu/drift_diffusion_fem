#include <iostream>
#include <string>

#include <armadillo>
#include "mesh.hpp"		// In "/home/xiaoxian/programs_cpp/triangle/

#include "exciton_diffusion.hpp"

using namespace std;

/************** Given a poly file:********************/
/* 1. Generate mesh on it by "Triangle"
/* 2. Assemble linear system 
/* 3. Solve the system
/*****************************************************/
int Solve(string polyname, const int num_squares, double& interface_drainage)
{
	/* Generate, plot, and save mesh */
	polyfileSmartInterface(polyname, num_squares);

	MeshData mesh;
	TriangleMesh(mesh, polyname);

	/* Linear system and solution */
	arma::mat M;
	arma::vec rhs;
	
	arma::vec a = interpolateCoefficient_a(mesh);
	arma::vec c = interpolateCoefficient_c(mesh);
	arma::vec d = interpolateCoefficient_d(mesh);
	arma::vec func_vec = interpolateRHS(mesh);

	arma::mat A = assemble_matrix_A(mesh, a);
	arma::mat C = assemble_matrix_C(mesh, c);
	arma::mat D = assemble_matrix_D(mesh, d);
	
	rhs = assemble_rhsvector(mesh, func_vec);
	M = A + C + D;
	

	applyDirichletBC(mesh, M, rhs);


	/* Solve and plot solution */
	arma::vec u;
	u = solve(M, rhs);

	PlotSolution(mesh, u, polyname);

/*	cout << "Matrix A:\n" << A << "\n";
	cout << "Matrix C:\n" << C << "\n";
	cout << "Matrix D:\n" << D << "\n";
	cout << "Matrix M:\n" << M << "\n";
*/

	interface_drainage = interfaceIntegral(mesh, u, d);
	cout << "Drainage on interface is " << interface_drainage << "\n";
	return 0;
}
