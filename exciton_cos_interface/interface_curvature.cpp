#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "../triangle/mesh.hpp"








double computeCurvature(const std::vector<double>& v0, const std::vector<double>& v1, const std::vector<double>& v2)
{
using namespace std;

	double curv = 0.0;

	double upper;
	double lower;

	upper = 0.5 * (v2[0] - v0[0]) * (v0[1] - 2*v1[1] + v2[1]) \
		-0.5 * (v0[0] - 2*v1[0] + v2[0]) * (v2[1] - v0[1]);
	lower = pow( 0.25 * ((v2[0] - v0[0]) * (v2[0] - v0[0]) + (v2[1] - v0[1]) * (v2[1] - v0[1])), 1.5);

	curv = upper / lower;
	return curv;
}


std::vector<double> interfaceCurvature(const my_mesh::MeshData &mesh)
{
using namespace std;
using namespace my_mesh;

	/* Check if interface edges have been generated */
	if (mesh.interface_edges.size()<=0)
	{	cout << "Need to compute interface nodes and edges first!\n";	exit(1);}

	/* Initialize curvature*/
	vector<double> curvatures;

	/* Compute curvature for each node on interface */
	for (int i=0; i<mesh.interface_nodes.size(); i++)
	{	vector<double> v0, v1, v2;
		if (i==0)
		{	int index1 = mesh.interface_nodes[i];
			int index2 = mesh.interface_nodes[i+1];

			v1 = mesh.nodes[index1];
			v2 = mesh.nodes[index2];

			v0.push_back( v2[0] );
			v0.push_back( 2.0 * v1[1] - v2[1] );
		}
		else if (i == mesh.interface_nodes.size()-1 )
		{	int index1 = mesh.interface_nodes[i];
			int index0 = mesh.interface_nodes[i-1];

			v1 = mesh.nodes[index1];
			v0 = mesh.nodes[index0];

			v2.push_back( v0[0] );
			v2.push_back( 2*v1[1] - v0[1] );
		}
		else
		{	int index0 = mesh.interface_nodes[i-1];
			int index1 = mesh.interface_nodes[i];
			int index2 = mesh.interface_nodes[i+1];

			v0 = mesh.nodes[index0];
			v1 = mesh.nodes[index1];
			v2 = mesh.nodes[index2];
		}

		int node_i_index = mesh.interface_nodes[i];
		double curv_i = computeCurvature(v0, v1, v2);
		curvatures.push_back(curv_i);
	}

	return curvatures;
}

