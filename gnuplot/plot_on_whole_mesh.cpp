#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

#include <armadillo>

#include "../triangle/mesh.hpp"

namespace my_gnuplot
{

/* Plot arma::vec over the entire domain "Omega" */
int plot_ArmaVec(const my_mesh::MeshData& mesh, const arma::vec u, const std::string& filename, std::string terminal)
{using namespace std;


	/* Write dat file for solution */
	ofstream output;
	string name_appended = "_whole_domain";
	string output_name = filename + name_appended + ".dat";
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
	string gnuplot_filename = filename + name_appended + ".gnuplot";
	ofstream gnuplot_fstream;
	gnuplot_fstream.open(gnuplot_filename.c_str());
	if (!gnuplot_fstream.is_open())
		{cout<< "Failed to open " << gnuplot_filename << "\n"; exit(1);	}

	if (terminal=="wxt")
	{	gnuplot_fstream << "set terminal wxt\n"
				<< "splot \"" << output_name << "\" with lines\n";
	}
	else if (terminal == "png")
	{	gnuplot_fstream	<< "set terminal png\n"
				<< "set output \"" << filename + name_appended + ".png" << "\"\n"
				<< "splot \"" << output_name << "\" with lines\n";
	}
	gnuplot_fstream.close();


	/* Run GNUplot file */
	string cmd = "gnuplot " + gnuplot_filename + " --persist\n";
	system(cmd.c_str());

	return 0;
}



}
