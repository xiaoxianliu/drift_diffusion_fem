#include <iostream>

#include "../mesh.hpp"

using namespace my_mesh;

int main(){
	std::string filename = "test_mesh";

	MeshData mesh_cos, mesh_greate_wall;
	mesh_cos = generateMesh_cosine_interface(filename, -0.5, 0.5, 2);
	mesh_greate_wall = generateMesh_great_wall(filename, -0.3, 1.0, 3);
	return 0;
}
