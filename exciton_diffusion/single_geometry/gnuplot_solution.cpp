#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <armadillo>

#include "mesh.hpp"
#include "exciton_diffusion.hpp"

using namespace std;

int PlotSolution(const MeshData& mesh, const arma::vec u, const string& polyname)
{	/* Write dat file for solution */
	ofstream output;
	string output_name = polyname + "_solution.dat";
	output.open(output_name.c_str());
	if (!output.is_open())
		{cout << "Failed to open file " << output_name << "\n"; exit(1);}

	/* Write coordinates and solution at each node; for loop through all elements */
	for (int t=0; t<mesh.num_elements; t++)
	{	int v0 = mesh.elements[3*t];
		int v1 = mesh.elements[3*t+1];
		int v2 = mesh.elements[3*t+2];

		output << mesh.nodes[2*v0] << "\t" << mesh.nodes[2*v0+1] << "\t" << u(v0) << "\n";
		output << mesh.nodes[2*v1] << "\t" << mesh.nodes[2*v1+1] << "\t" << u(v1) << "\n";
		output << mesh.nodes[2*v2] << "\t" << mesh.nodes[2*v2+1] << "\t" << u(v2) << "\n";
		output << mesh.nodes[2*v0] << "\t" << mesh.nodes[2*v0+1] << "\t" << u(v0) << "\n\n\n";
	}
	output.close();

	/* Write GNUplot file */
	string gnuplot_filename = polyname + "_solution.gnuplot";
	ofstream gnuplot_fstream;
	gnuplot_fstream.open(gnuplot_filename.c_str());
	if (!gnuplot_fstream.is_open())
		{cout<< "Failed to open " << gnuplot_filename << "\n"; exit(1);	}

	gnuplot_fstream	<< "set terminal png\n"
			<< "set output \"" << polyname << "_solution.png\"\n"
			<< "splot \"" << output_name << "\" with lines\n";

//	gnuplot_fstream << "splot \"" << output_name << "\" with lines\n";

	gnuplot_fstream.close();


	/* Run GNUplot file */
//	string cmd = "gnuplot " + gnuplot_filename + "\n";
	string cmd = "gnuplot " + gnuplot_filename + " --persist\n";
	system(cmd.c_str());

	return 0;
}
