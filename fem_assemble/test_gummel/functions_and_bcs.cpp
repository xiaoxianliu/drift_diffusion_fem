#include <vector>
#include <armadillo>
#include "../../triangle/mesh.hpp"

double psi_func(double x, double y)
{	return 1+x+y;
}

double f_func(double x, double y)
{	return -1-x;
}






// Boundary conditions
double uD_func(double x, double y)
{	return 0.5;
}


int applyDirichletBC(	const my_mesh::MeshData &mesh,
			arma::mat &M,
			arma::vec &rhs)
{

	for (int i=0; i<mesh.num_nodes; i++)
	{	std::vector<int> neigh_edges = mesh.topology0to1[i];

		bool is_dirichlet = false;
		for (int e=0; e<neigh_edges.size(); e++)
		{	if ( mesh.edge_markers[ neigh_edges[e] ] == 1 || mesh.edge_markers[ neigh_edges[e] ] == 3 )
			{	is_dirichlet = true;
				break;
			}
		}

		if (is_dirichlet)
		{	double x = mesh.nodes[i][0], y = mesh.nodes[i][1];
			M(i, arma::span::all).zeros();
			M(i,i) = 1.0;
			rhs(i) = uD_func(x,y);
		}
	}
	return 0;
}









// plot solutions
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
