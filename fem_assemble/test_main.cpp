#include <vector>
#include <armadillo>

#include "fem_assemble.hpp"
#include "../triangle/mesh.hpp"




my_mesh::MeshData generateMesh();
double rhsFunc(double x, double y);
arma::vec rhsVec(const my_mesh::MeshData& mesh);
double uD_left(double x, double y);
double uD_right(double x, double y);
int applyDirichletBC(const my_mesh::MeshData& mesh, arma::mat& M, arma::vec& b);


int main()
{
	using namespace my_mesh;
//	using namespace linear_fem;

	/* 1. mesh */
	MeshData mesh;
	mesh = generateMesh();

	/* 2. coefficient matrix */
	arma::vec a = arma::ones<arma::vec>(mesh.num_nodes);
	arma::mat A = assembleMatrixA(mesh,a);

	return 0;
}




my_mesh::MeshData generateMesh()
{
using namespace std;
using namespace my_mesh;

	/* 1. Form new mesh input */
	vector<vector<double> > interface_nodes;
	vector<int> node_markers;
	vector<vector<int> > segments;
	vector<int> segment_markers;
	vector<vector<double> > regions;
	vector<int> region_markers;

	/* 1.1 Add interface node to the vector of interface nodes */
	vector<double> new_node(2);
	new_node[0] = -0.5;	new_node[1]=1.0;
	interface_nodes.push_back(new_node);
	new_node[0] = 0.2;	new_node[1]=0.5;
	interface_nodes.push_back(new_node);
	new_node[0] = -0.3;	new_node[1] = 0.0;
	interface_nodes.push_back(new_node);


	newMesh(interface_nodes, node_markers,		// 4 corner vertices are added to "interface_nodes"
		segments, segment_markers,
		regions, region_markers);

	vector<vector<double> > nodes = interface_nodes;

	/* 2. Write to poly file */
	string filename="test_fem";
	writePolyfile(	filename,
			nodes, node_markers,
			segments, segment_markers,
			regions, region_markers);


	/* 3. Generate mesh by "triangle" */
	string cmd = "/home/xiaoxian/bin/triangle/triangle -qzpeAa.002 " + filename + ".poly";
	system(cmd.c_str());

	/* 4. read in mesh information */
	MeshData mesh;

	string node_file_name= filename + ".1.node";
	string edge_file_name= filename + ".1.edge";
	string ele_file_name = filename + ".1.ele";

	ReadNodes(mesh, node_file_name);
	ReadEdges(mesh, edge_file_name);
	ReadElements(mesh, ele_file_name);

	ComputeTopology(mesh);

	ComputeEdgeLengths(mesh);
	ComputeElementAreas(mesh);

	return mesh;	
}

















// Right-hand side function
double rhsFunc(double x, double y)
{	return -1.0+0.5*x*x;
}

// Right-hand side vector (projected to linear FEM space)
arma::vec rhsVec(const my_mesh::MeshData& mesh)
{
using namespace arma;
using std::vector;
	int num_nodes = mesh.num_nodes;

	vec b(num_nodes);
	b.zeros();

	for (int v; v<num_nodes; v++)
	{
		double x = mesh.nodes[v][0];
		double y = mesh.nodes[v][1];

		vector<int> Ts = mesh.topology0to2[v];
		for (int i=0; i<Ts.size(); i++)
		{	int t = Ts[i];
			double area_T = mesh.ele_areas[t];
			b(v) += rhsFunc(x,y) * area_T / 3.0;
		}

	}

}

// Functions defining Dirichlet functions
double uD_left(double x, double y)
{	return 0.5;	}
double uD_right(double x, double y)
{	return 0.5;	}


/* Impose Dirichlet boundary conditions */
/* By default, "1" and "3" are the "left" and "right" boundaries denoting Dirichlet boundary conditions */
int applyDirichletBC(const my_mesh::MeshData& mesh, arma::mat& M, arma::vec& b)
{
using std::vector;
using namespace arma;

	int num_nodes = mesh.num_nodes;
	for (int v; v<num_nodes; v++)
	{
		vector<int> Es = mesh.topology0to1[v];
		for (int i=0; i<Es.size(); i++)
		{	int edge = Es[i];
			int edge_marker = mesh.edge_markers[edge];
			if (edge_marker == 1)			// left boundary condition
			{	double x = mesh.nodes[v][0];
				double y = mesh.nodes[v][1];

				M(v, span::all) = zeros<mat>(1, num_nodes);
				M(v,v) = 1.0;
				b(v) = uD_left(x,y);
			}
			else if (edge_marker == 3)
			{	double x = mesh.nodes[v][0];
				double y = mesh.nodes[v][1];

				M(v, span::all) = zeros<mat>(1, num_nodes);
				M(v,v) = 1.0;
				b(v) = uD_right(x,y);
			}
		}
	}
}
