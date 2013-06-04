#include <vector>
#include <string>
#include <armadillo>

#include "../fem_assemble.hpp"
#include "../../triangle/mesh.hpp"




my_mesh::MeshData generateMesh(std::string filename);
double rhsFunc(double x, double y);
arma::vec rhsVec(const my_mesh::MeshData& mesh);
double uD_left(double x, double y);
double uD_right(double x, double y);
int applyDirichletBC(const my_mesh::MeshData& mesh, arma::mat& M, arma::vec& b);
int plotSolution(const my_mesh::MeshData& , const arma::vec u, const std::string& filename);


int main()
{
//	using namespace std;
	using namespace my_mesh;
	using namespace linear_fem;


	std::string filename = "test_fem";

	/* 1. mesh */
	MeshData mesh;
	mesh = generateMesh(filename);

	/* 2. coefficient matrix */

	arma::mat A = assembleMatrixA(mesh, arma::ones<arma::vec>(mesh.num_nodes));
	arma::mat C = assembleMatrixC(mesh, arma::ones<arma::vec>(mesh.num_nodes));
	arma::mat M = A + C;
	arma::vec b = rhsVec(mesh);
	applyDirichletBC(mesh, M, b);

	/* 3. Solve and plot solutions */
	arma::vec u;
	u = arma::solve(M,b);
	plotSolution(mesh, u, filename);

	return 0;
}











/************************* Function related to this test ***************************/

/* Calling "triangle" to generate mesh */

my_mesh::MeshData generateMesh(std::string filename)
{
using namespace my_mesh;
using namespace std;

	/* 1. Form new mesh input */
	vector<vector<double> > interface_nodes;
	vector<int> node_markers;
	vector<vector<int> > segments;
	vector<int> segment_markers;
	vector<vector<double> > regions;
	vector<int> region_markers;

	/* 1.1 Add interface node to the vector of interface nodes */
	vector<double> new_node(2);
	new_node[0] = 0.0;	new_node[1]=1.0;
	interface_nodes.push_back(new_node);
	new_node[0] = 0.0;	new_node[1]=0.0;
	interface_nodes.push_back(new_node);


	newMesh(interface_nodes, node_markers,		// 4 corner vertices are added to "interface_nodes" inside "newMesh(...)"
		segments, segment_markers,
		regions, region_markers);

	vector<vector<double> > nodes = interface_nodes;

	/* 2. Write to poly file */
//	string filename="test_fem";
	writePolyfile(	filename,
			nodes, node_markers,
			segments, segment_markers,
			regions, region_markers);


	/* 3. Generate mesh by "triangle" */
	string cmd = "/home/xiaoxian/bin/triangle/triangle -qzpeAa0.005 " + filename + ".poly";
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

	gnuplot_mesh(mesh, filename);
	return mesh;	
}









/*  Right-hand side function  */
double rhsFunc(double x, double y)
{	return -1.0+0.5*x*x;
}

/* Right-hand side vector after projection */
arma::vec rhsVec(const my_mesh::MeshData& mesh)
{
using namespace arma;
using std::vector;
	int num_nodes = mesh.num_nodes;

	vec b(num_nodes);
	b.zeros();

	for (int v=0; v<num_nodes; v++)
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
	return b;

}






/* Functions defining Dirichlet functions */
double uD_left(double x, double y)
{	return 0.5;	}
double uD_right(double x, double y)
{	return 0.5;	}






/* Impose Dirichlet boundary conditions */
/* By default, "1" and "3" are the "left" and "right" boundaries denoting Dirichlet boundary conditions */
int applyDirichletBC(const my_mesh::MeshData& mesh, arma::mat& M, arma::vec& b)
{

using namespace std;
using namespace arma;

	int num_nodes = mesh.num_nodes;
	for (int v=0; v<num_nodes; v++)
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
			else if (edge_marker == 3)		// right boundary condition
			{	double x = mesh.nodes[v][0];
				double y = mesh.nodes[v][1];

				M(v, span::all) = zeros<mat>(1, num_nodes);
				M(v,v) = 1.0;
				b(v) = uD_right(x,y);
			}
		}

	}

	return 0;
}






/* Plot solution */
int plotSolution(const my_mesh::MeshData& mesh, const arma::vec u, const std::string& filename)
{using namespace std;


	/* Write dat file for solution */
	ofstream output;
	string output_name = filename + "_solution.dat";
	output.open(output_name.c_str());
	if (!output.is_open())
		{cout << "Failed to open file " << output_name << "\n"; exit(1);}

	/* Write coordinates and solution at each node; for loop through all elements */
	for (int t=0; t<mesh.num_elements; t++)
	{	int v0 = mesh.elements[t][0];
		int v1 = mesh.elements[t][1];
		int v2 = mesh.elements[t][2];

		output << mesh.nodes[v0][0] << "\t" << mesh.nodes[v0][1] << "\t" << u(v0) << "\n";
		output << mesh.nodes[v1][0] << "\t" << mesh.nodes[v1][1] << "\t" << u(v1) << "\n";
		output << mesh.nodes[v2][0] << "\t" << mesh.nodes[v2][1] << "\t" << u(v2) << "\n";
		output << mesh.nodes[v0][0] << "\t" << mesh.nodes[v0][1] << "\t" << u(v0) << "\n\n\n";
	}
	output.close();

	/* Write GNUplot file */
	string gnuplot_filename = filename + "_solution.gnuplot";
	ofstream gnuplot_fstream;
	gnuplot_fstream.open(gnuplot_filename.c_str());
	if (!gnuplot_fstream.is_open())
		{cout<< "Failed to open " << gnuplot_filename << "\n"; exit(1);	}

/*	gnuplot_fstream	<< "set terminal png\n"
			<< "set output \"" << filename << "_solution.png\"\n"
			<< "splot \"" << output_name << "\" with lines\n";
*/
	gnuplot_fstream << "splot \"" << output_name << "\" with lines\n";

	gnuplot_fstream.close();


	/* Run GNUplot file */
//	string cmd = "gnuplot " + gnuplot_filename + "\n";
	string cmd = "gnuplot " + gnuplot_filename + " --persist\n";
	system(cmd.c_str());

	return 0;
}
