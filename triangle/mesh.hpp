/* This file construct a class/struct called "mesh" to store various information of the mesh generated by "triangle"
*/
#include <vector>
#include <string>



#ifndef MESH_H
#define MESH_H

#define MeshNamespace my_mesh				// "my_mesh" is replaced with "MeshNamespace" throughout this folder

namespace MeshNamespace
{

/* Structure for storing mesh information (e.g. topology, triangle areas, etc.)*/
struct MeshData;

/********* Function declarations ****************/
/* New mesh input */
int newMesh(	std::vector< std::vector<double> >& nodes, std::vector<int>& node_markers,\
		std::vector<std::vector<int> >& segments,	std::vector<int>& segment_markers,\
		std::vector<std::vector<double> >& regions, std::vector<int>& region_markers);


/* Read in mesh information and write to a .poly file as input for "triangle" */
int writePolyfile(std::string filename,
		std::vector< std::vector<double> > nodes, std::vector<int> node_markers,		// vertices
		std::vector< std::vector<int> > segments, std::vector<int> segment_markers,		// segments
		std::vector< std::vector<double> > regions, std::vector<int> region_markers);		// regional attributes


/* Read mesh data from ".node", ".edge", ".ele" files */
int ReadNodes(MeshData&, std::string);
int ReadEdges(MeshData&, std::string);
int ReadElements(MeshData&, std::string);

/* vector operations */
int add_one_entry(std::vector<int>&, const int&);
int subtract_one_entry(std::vector<int>&, const int&);
std::vector<int> intersect_vectors(const std::vector<int>&, const std::vector<int>&);
std::vector<int> merge_vectors(const std::vector<int>&, const std::vector<int>&);

/* Compute mesh topology (i.e. connectivity) */
int ComputeTopology(MeshData &mesh);

/* Compute geometric properties */
void ComputeElementAreas(MeshData &mesh);
void ComputeEdgeLengths(MeshData &mesh);
void ComputeInterfaceCurvatures(MeshData &mesh);

/* Write .gnuplot file to plot mesh */
void gnuplot_mesh(const MeshData &mesh, const std::string& filename);
void gnuplot_interface(const MeshData &mesh, const std::string& filename);

/* Find sub-mesh or interface */
int extractInterface(	const MeshData& mesh, std::vector<int>& interface_edges, std::vector<int>& interface_nodes);

}





// Define struct "MeshData"
struct MeshNamespace::MeshData
{
	/* Nodes */
	size_t num_nodes, num_attributes_per_node, num_marker_per_node;
	std::vector< std::vector<double> > nodes;			// coordinates of nodes (num_nodes * 2)
	std::vector< std::vector<double> > node_attributes;	// attributes of each node; each node could have multiple attributes.
	std::vector<int> node_markers;

	/* Edges */
	size_t num_edges;
	size_t num_marker_per_edge;
	std::vector< std::vector<int> > edges;			// indices of end nodes of edges (num_edges * 2)
	std::vector<int> edge_markers;

	/* Elements */
	size_t num_elements;
	size_t num_nodes_per_ele;
	size_t num_attributes_per_ele;
	std::vector< std::vector<int> > elements;		// indices of all nodes of element (num_elements * num_nodes_per_ele)
	std::vector<int> element_markers;			// Assume this is the only attribute each element has!!!

	/* Interface and its (nodal) curvatures*/
	std::vector< int > interface_nodes;			// indices of nodes and edges on interface
	std::vector< int > interface_edges;			// Note: num_interface_nodes = num_interface_edges+1
	std::vector<double> interface_curvatures;

	/* Mesh topology */
	std::vector< std::vector<int> > topology2to0;	// nodes of each element (3 for linear element, 6 for quadratic element)
	std::vector< std::vector<int> > topology2to1;		// edges of each element (always 3)
	std::vector< std::vector<int> > topology2to2;		// neighboring element (always 3)

	std::vector< std::vector<int> > topology1to0;		//nodes of edge (always 2)
	std::vector< std::vector<int> > topology1to1;		//edges sharing nodes 
	std::vector< std::vector<int> > topology1to2;		//elements containing edge (always 2)

	std::vector< std::vector<int> > topology0to0;		//neighboring nodes
	std::vector< std::vector<int> > topology0to1;		//edges containing node
	std::vector< std::vector<int> > topology0to2;		//elements containing node

	/* Mesh properties */
	std::vector<double> edge_lengths;
	std::vector<double> ele_areas;





	int clear_meshdata()
	{
		num_nodes=0; num_attributes_per_node=0;	num_marker_per_node=0;
		nodes.clear();	
		node_attributes.clear();
		node_markers.clear();

		num_edges=0;	num_marker_per_edge=0;
		edges.clear();	
		edge_markers.clear();

		/* Elements */
		num_elements=0;	
		num_nodes_per_ele=0;
		num_attributes_per_ele=0;
		elements.clear();	
		element_markers.clear();

		/* Interface */
		interface_nodes.clear();
		interface_edges.clear();

		/* Mesh topology */
		topology2to0.clear();	
		topology2to1.clear();	
		topology2to2.clear();	

		topology1to0.clear();	
		topology1to1.clear();	
		topology1to2.clear();	

		topology0to0.clear();	
		topology0to1.clear();	
		topology0to2.clear();	

		/* Mesh properties */
		edge_lengths.clear();
		ele_areas.clear();

		return 0;
	}
	
};







#endif
