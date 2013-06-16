#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

#include <armadillo>

#include "../triangle/mesh.hpp"


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
