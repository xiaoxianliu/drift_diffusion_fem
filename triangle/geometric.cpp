#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "mesh.hpp"

using namespace std;
using namespace MeshNamespace;

/***** 1. Compute element area *****/

void MeshNamespace::ComputeElementAreas(MeshData &mesh){
	mesh.ele_areas.clear();
	mesh.ele_areas.resize(mesh.num_elements);

	for (int i=0; i<mesh.num_elements; i++)
	{	/* Node indices of element "i"*/
		int v0 = mesh.elements[i][0];
		int v1 = mesh.elements[i][1];
		int v2 = mesh.elements[i][2];

		/* Compute 2 edges of element "i"*/
		double edge0[2], edge1[2];
		for (int j=0; j<2; j++)
		{	edge0[j] = mesh.nodes[v1][j] - mesh.nodes[v0][j];
			edge1[j] = mesh.nodes[v2][j] - mesh.nodes[v0][j];
		}

		/* Compute area_i (=|edge0 x edge1|/2) */
		mesh.ele_areas[i] = 0.5 * (edge0[0]*edge1[1] - edge0[1]*edge1[0]);
	}
}


/***** 2. Compute edge length *****/
void MeshNamespace::ComputeEdgeLengths(MeshData &mesh){
	mesh.edge_lengths.clear();
	mesh.edge_lengths.resize(mesh.num_edges);

	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[i][0];
		int v1 = mesh.edges[i][1];

		double edge[2];
		edge[0]= mesh.nodes[v1][0] - mesh.nodes[v0][0];
		edge[1]= mesh.nodes[v1][1] - mesh.nodes[v0][1];

		mesh.edge_lengths[i] = sqrt(edge[0]*edge[0] + edge[1]*edge[1]);
	}
}




/***** 3. Compute interface curvatures *********/

// Local function: estimate the curvature at given node of some curve 
double computeCurvature(const std::vector<double>& v0, 
			const std::vector<double>& v1, 
			const std::vector<double>& v2)
{
	double curv = 0.0;

	double upper;
	double lower;

	upper = 0.5 * (v2[0] - v0[0]) * (v0[1] - 2*v1[1] + v2[1]) \
		-0.5 * (v0[0] - 2*v1[0] + v2[0]) * (v2[1] - v0[1]);
	lower = pow( 0.25 * ((v2[0] - v0[0]) * (v2[0] - v0[0]) + (v2[1] - v0[1]) * (v2[1] - v0[1])), 1.5);

	curv = upper / lower;
	return curv;
}




/*** Compute curvature on the entire interface ***/
void MeshNamespace::ComputeInterfaceCurvatures( MeshData &mesh )
{
	/* Check if interface edges have been generated */
	if (mesh.interface_edges.size()<=0)
	{	cout << "Need to compute interface nodes and edges first!\n";	exit(1);}

	mesh.interface_curvatures.clear();

	/* Compute curvature for each node on interface */
	for (int i=0; i<mesh.interface_nodes.size(); i++)
	{	vector<double> v0, v1, v2;
		if (i==0)					//i=0 corresponds to first node (top boundary)
		{	int index1 = mesh.interface_nodes[i];
			int index2 = mesh.interface_nodes[i+1];

			v1 = mesh.nodes[index1];
			v2 = mesh.nodes[index2];

			v0.push_back( v2[0] );
			v0.push_back( 2.0 * v1[1] - v2[1] );
		}
		else if (i == mesh.interface_nodes.size()-1 )	//i=interface_nodes.size()-1 corresponds to last node (bottome bdry)
		{	int index1 = mesh.interface_nodes[i];
			int index0 = mesh.interface_nodes[i-1];

			v1 = mesh.nodes[index1];
			v0 = mesh.nodes[index0];

			v2.push_back( v0[0] );
			v2.push_back( 2*v1[1] - v0[1] );
		}
		else						// nodes in the middle
		{	int index0 = mesh.interface_nodes[i-1];
			int index1 = mesh.interface_nodes[i];
			int index2 = mesh.interface_nodes[i+1];

			v0 = mesh.nodes[index0];
			v1 = mesh.nodes[index1];
			v2 = mesh.nodes[index2];
		}

		double curv_i = computeCurvature(v0, v1, v2);
		mesh.interface_curvatures.push_back(curv_i);
	}
}

