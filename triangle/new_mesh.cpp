#include <iostream>
#include <cstdlib>
#include <vector>
#include "mesh.hpp"

#define EOF 5e-16

using namespace std;
using namespace MeshNamespace;



/* Given information of interface nodes in "nodes" as input, form a PSLG plot after adding (-1,0), (1,0), (1,1), (-1,1) to the list of nodes */
/* Nonempty input: only interface "nodes" */
int MeshNamespace::initializeMesh(	vector< vector<double> >& nodes,\
					vector<int>& node_markers,\
					vector<vector<int> >& segments,\
					vector<int>& segment_markers,\
					vector<vector<double> >& regions,\
					vector<int>& region_markers)
{
	/* 1. Form list of nodes */
	int num_interface_nodes = nodes.size();
	if (num_interface_nodes < 2)
		{cout<< "Num of interface nodes is " << nodes.size()<<"; it should be at least 2! Exit...\n";	exit(1);}

	/* 1.1 Add corner node coordinates to the vector of "nodes" */
	vector<vector<double> > corner_nodes(4);
	corner_nodes[0].push_back(1);	corner_nodes[0].push_back(0);
	corner_nodes[1].push_back(1);	corner_nodes[1].push_back(1);
	corner_nodes[2].push_back(-1);	corner_nodes[2].push_back(1);
	corner_nodes[3].push_back(-1);	corner_nodes[3].push_back(0);
	nodes.insert(nodes.begin(), corner_nodes.begin(), corner_nodes.end());

	/* 1.2 Add corner node markers if the given "node_markers" is not empty */
	if (node_markers.size() == num_interface_nodes)
	{	vector<int> corner_node_markers(4);
		for (int i=0; i<corner_node_markers.size(); i++)
		{	corner_node_markers[i] = 0;					// set markers of corner nodes to "0"
		}
		node_markers.insert( node_markers.begin(), corner_node_markers.begin(), corner_node_markers.end());
	}

	/* 2. Form the segment */
	segments.clear();	segment_markers.clear();

	/* 2.1 Identify which end of the interface is on the top boundary */
	int top_end = -1, bottom_end = -1;
	vector<double> endnode1 = nodes[4], endnode2 = nodes[nodes.size()-1];
	if (endnode1[1]>1-EOF && endnode2[1]<EOF)
		{top_end = 4;	bottom_end = nodes.size()-1;	}
	else if (endnode1[1]<EOF && endnode2[1]>1-EOF)
		{top_end = nodes.size()-1;	bottom_end = 4;	}
	else
	{  cout << "Not both ends of interface are on the boundary of the 2-by-1 rectangle!"; exit(1);}

	/* 2.2 Form the segments on outer boundaries */
	vector<int> boundary_node_indices(7);
	vector<int> boundary_segment_markers(6);

						/* Boundary markers: left 1, bottom 2, right 3, top 4. */
	boundary_node_indices[0] = 0;
	boundary_node_indices[1] = 1;		boundary_segment_markers[0] = 3;
	boundary_node_indices[2] = top_end;	boundary_segment_markers[1] = 4;
	boundary_node_indices[3] = 2;		boundary_segment_markers[2] = 4;
	boundary_node_indices[4] = 3;		boundary_segment_markers[3] = 1;
	boundary_node_indices[5] = bottom_end;	boundary_segment_markers[4] = 2;
	boundary_node_indices[6] = 0;		boundary_segment_markers[5] = 2;

	vector<int> new_seg(2);
	for (int i=0; i<6; i++)
	{
		new_seg[0] = boundary_node_indices[i];	new_seg[1] = boundary_node_indices[i+1];
		segments.push_back(new_seg);
		segment_markers.push_back (boundary_segment_markers[i]);
	}

	/* 2.3 Form the segments on interface */
	for (int i=4; i<nodes.size()-1; i++)
	{	new_seg[0] = i;	new_seg[1] = i+1;
		segments.push_back(new_seg);
		segment_markers.push_back(5);						// interface segments are marked by "5"
	}


	/* 3. Form regions (no input) */
	regions.clear();	region_markers.clear();

	vector<double> new_reg(2);
	new_reg[0] = -0.99;	new_reg[1] = 0.99;					// Left region
	regions.push_back(new_reg);	region_markers.push_back(1);
	new_reg[0] = 0.99;	new_reg[1] = 0.99;
	regions.push_back(new_reg);	region_markers.push_back(2);			// Right region

	return 0;

}
