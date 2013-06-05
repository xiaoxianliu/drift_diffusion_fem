#include <string>

#include "exciton_cos_interface.hpp"
#include "../triangle/mesh.hpp"

int main()
{
	std::string filename = "cos_interface";
	double y_control = 0.5;

	/* 1. Mesh */
	my_mesh::MeshData mesh;
	mesh = generateMesh_CosInterface(filename, y_control);

	return 0;
}
