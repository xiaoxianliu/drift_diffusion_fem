#include <iostream>
#include <vector>
#include <string>
#include <armadillo>

#include "../fem_assemble.hpp"
#include "../../triangle/mesh.hpp"
#include "test.hpp"


#define PI 3.1415926535897932384626433832795



int main()
{
using namespace std;
using namespace my_mesh;

	std::string filename = "convergence";

	/* 0. Max areas for a sequence of meshes */
	std::vector<double> max_areas;
	max_areas.push_back(0.1);
	max_areas.push_back(0.03);
	max_areas.push_back(0.01);
	max_areas.push_back(0.003);
	max_areas.push_back(0.001);


	/* 1.1 Inteface node coordinates */
	int num_interface_nodes = 501;
	vector< vector<double> > interface_nodes;

	double y_control = 0.5;
	for (int i=0; i<num_interface_nodes; i++)
	{	vector<double> node(2);
		double y;
		if (i<num_interface_nodes-1)
			y = 1.0 - i/static_cast<double>(num_interface_nodes-1);
		else
			y = 0.0;
		double x = y_control * 0.5*(1 - std::cos(2*PI*y));

		node[0] = x;
		node[1] = y;

		interface_nodes.push_back(node);
	}


	/* 2. Generate sequence of meshes and then solve for solution*/
	for (int i=0; i<max_areas.size(); i++)
	{
		std::cout << "max_area[" << i << "] = " << max_areas[i] << "\n";

		/* 1. mesh */
		MeshData mesh;
		mesh = generateMesh(filename, interface_nodes, max_areas[i]);

		/* 2. Solve for solution and error */
		arma::vec u(mesh.num_nodes);
		double error;
		solveEq(mesh, filename, u, error);

		std::cout << "maximum error from true solution is " << error << "\n";
	}






//	double interface_loss = interfaceIntegral(mesh,u);
//	std::cout<< "interface_loss is " << interface_loss << "\n";

	return 0;
}










