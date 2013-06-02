#include <iostream>
#include <cstdlib>
#include <vector>
#include "mesh.hpp"
using namespace std;
using namespace MeshNamespace;

#define EOF 5.0e-16

int MeshNamespace::extractInterface(	const MeshData& mesh, \
					vector<int>& interface_edges, \
					vector<int>& interface_nodes)
{
	interface_edges.clear();
	interface_nodes.clear();

	/* Find first node on "top" segment: the node = (x,y) with y==1*/
	int first_node=-1;
	for (int i=0; i<mesh.num_edges; i++)
	{	if (mesh.edge_markers[i]==5)				// By assumption, interface segments are marked with "5"
		{	int v0 = mesh.edges[i][0];
			int v1 = mesh.edges[i][1];
			double y0 = mesh.nodes[v0][1];
			double y1 = mesh.nodes[v1][1];
			if ( (y0-1.0)<EOF && (y0-1.0)>-EOF ) 
				first_node = v0;
			if ( (y1-1.0)<EOF && (y1-1.0)>-EOF )
				first_node = v1;
		}
		if (first_node >=0)	break;
	}

	if (first_node >=0)
	{	interface_nodes.push_back(first_node);
	}
	else
	{	cout << "Didn't find the node on top segment with y-coordinate == 1... Go check where went wrong!\n"; exit(1);}


	/* Next, recover the whole "chain" of edges and nodes following first_node */
	int last_node = first_node;					// index of the last node of interface chain
	int new_node;
	double last_node_y = mesh.nodes[last_node][1];			// y-coordinate of the last node of interface chain
	int last_edge = -1;						// index of the last edge of inteface chain
	int new_edge = -1;
	bool new_edge_found = true;

	while (last_node_y > EOF && new_edge_found)
	{	new_edge_found = false;

		/* Search for new neighboring edge of "last_node" , and add that to "interface_edge"*/
		vector<int> neigh_edges = mesh.topology0to1[last_node];
		for (int i=0; i<neigh_edges.size(); i++)
		{	int edge = neigh_edges[i];
			if (mesh.edge_markers[edge] == 5 && edge != last_edge)			// if a new neighboring edge is found
			{	new_edge = edge;
				last_edge = new_edge;
				interface_edges.push_back(new_edge);

				new_edge_found = true;
				break;
			}
		}

		/* Add newly found edge and nodes to "interface_edges" and "interface_nodes" */
		if (new_edge_found)
		{	int v0 = mesh.edges[new_edge][0];
			int v1 = mesh.edges[new_edge][1];
			if (last_node == v0)	new_node = v1;
			else if (last_node == v1)	new_node = v0;
			else
			{cout<< "Nodes of edge " << new_edge << " are " << v0 << " " << v1 << "\n";
			 cout<< "last node of existing interface_chain is " << last_node << "\n";
			 cout<< "---> Last node must be one of " << v0 << " and " << v1 << "\n";
			}

			interface_nodes.push_back(new_node);
			last_node = new_node;
			last_node_y = mesh.nodes[new_node][1];
		}
	}

	return 0;
}





