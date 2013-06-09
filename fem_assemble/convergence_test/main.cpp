#include <iostream>
#include <vector>
#include <string>
#include <armadillo>

#include "../fem_assemble.hpp"
#include "../../triangle/mesh.hpp"
#include "test.hpp"






int main()
{

	using namespace my_mesh;
	std::string filename = "convergence";

	/* 0. Max areas for a sequence of meshes */
	std::vector<double> max_areas;
	max_areas.push_back(0.1);
	max_areas.push_back(0.03);
	max_areas.push_back(0.01);
	max_areas.push_back(0.003);
	max_areas.push_back(0.001);

	for (int i=0; i<max_areas.size(); i++)
	{
		std::cout << "max_area[" << i << "] = " << max_areas[i] << "\n";

		/* 1. mesh */
		MeshData mesh;
		mesh = generateMesh_CosInterface(filename, -0.5, max_areas[i]);

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










