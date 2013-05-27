#include "mesh.hpp"
#include "cmath"
using namespace std;

/* Compute element area and edge length */

void ComputeElementAreas(MeshData &mesh){
	mesh.ele_areas = new double[mesh.num_elements];

	for (int i=0; i<mesh.num_elements; i++)
	{	/* Node indices of element "i"*/
		int v0 = mesh.elements[3*i];
		int v1 = mesh.elements[3*i+1];
		int v2 = mesh.elements[3*i+2];

		/* Compute 2 edges of element "i"*/
		double edge0[2], edge1[2];
		for (int j=0; j<2; j++)
		{	edge0[j] = mesh.nodes[2*v1+j] - mesh.nodes[2*v0+j];
			edge1[j] = mesh.nodes[2*v2+j] - mesh.nodes[2*v0+j];
		}

		/* Compute area_i (=|edge0 x edge1|/2) */
		mesh.ele_areas[i] = 0.5 * (edge0[0]*edge1[1] - edge0[1]*edge1[0]);
	}
}

void ComputeEdgeLengths(MeshData &mesh){
	mesh.edge_lengths = new double[mesh.num_edges];

	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[2*i];
		int v1 = mesh.edges[2*i+1];

		double edge[2];
		edge[0]= mesh.nodes[2*v1] - mesh.nodes[2*v0];
		edge[1]= mesh.nodes[2*v1+1] - mesh.nodes[2*v0+1];

		mesh.edge_lengths[i] = sqrt(edge[0]*edge[0] + edge[1]*edge[1]);
	}
}
