#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include <armadillo>

#include "exciton_interface_flux.hpp"
#include "../triangle/mesh.hpp"

#define PI 3.1415926535897932384626433832795

int main()
{
using namespace std;

	std::string filename = "Cosine";

	/* 1. Mesh */
	/* 1.1 Generate sample points on interface */
	// all sample points with coordinates (x,y) satisfy the formula:
	// x = y_control * (1 - cos(2*pi*y)); 
	int num_interface_nodes = 401;
	vector< vector<double> > interface_nodes;

	double x_control = 0.5;
	for (int i=0; i<num_interface_nodes; i++)
	{	vector<double> node(2);
		double y;
		if (i<num_interface_nodes-1)
			y = 1.0 - i/static_cast<double>(num_interface_nodes-1);
		else
			y = 0.0;
		double x = x_control * 0.5*(1 - std::cos(2*PI*y));

		node[0] = x;
		node[1] = y;

		interface_nodes.push_back(node);
	}

	/* 1.2 Generate mesh */
	my_mesh::MeshData mesh;
	double max_area = 0.001;
	mesh = my_mesh::generateMesh(filename, interface_nodes, max_area);

	/* 1.3 Compute curvature of interface and plot it */
	my_mesh::ComputeInterfaceCurvatures(mesh);
	plot_InterfaceSTLVector(mesh, mesh.interface_curvatures, "curvatures");	// plot interface curvature




	// 2. State equation 
	arma::vec u = solveStateEq(mesh);
	plot_ArmaVec(mesh, u, filename+"_state");
	plot_ArmaVec_on_Interface(mesh, u, filename+"_state");			// plot State Eq solution on interface


	// 3. Compute normal derivative on interface
	arma::vec du_dnu_1, du_dnu_2;
	computeInterfaceNormalDerivative (mesh, u, du_dnu_1, du_dnu_2);
	plot_InterfaceArmaVec(mesh, du_dnu_1, filename+"_interface_normal_deriv_1");
	plot_InterfaceArmaVec(mesh, du_dnu_2, filename+"_interface_normal_deriv_2");



	/* 5. Compute functional */
	double loss = interfaceIntegral( mesh, u );
	std::cout << "total loss at interface is " << loss << "\n";





	return 0;
}
