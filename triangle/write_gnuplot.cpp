#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include "mesh.hpp"
using namespace std;


void WriteGNUplot(MeshData &mesh, const string& polyname){

	/* Write data file */
	ofstream output;
	string filename = polyname + "_mesh.dat";
	output.open(filename.c_str());

	if (!output.is_open())
		{cout << "Fail to open "<<filename<<endl; exit(1);}

	for (int i=0; i<mesh.num_elements; i++)
	{	int v0 = mesh.elements[3*i];
		int v1 = mesh.elements[3*i+1];
		int v2 = mesh.elements[3*i+2];

		output<<mesh.nodes[2*v0] << "\t" << mesh.nodes[2*v0+1] << "\n";
		output<<mesh.nodes[2*v1] << "\t" << mesh.nodes[2*v1+1] << "\n";
		output<<mesh.nodes[2*v2] << "\t" << mesh.nodes[2*v2+1] << "\n";
		output<<mesh.nodes[2*v0] << "\t" << mesh.nodes[2*v0+1] << "\n\n";
	}

	output.close();

	/* Write GNUplot command */
	string meshfile_name;
	meshfile_name = polyname + "_mesh.gnuplot";
	ofstream gnuplot_file;
	gnuplot_file.open(meshfile_name.c_str());
	if (!gnuplot_file.is_open())								// exit if failed to open file
		{cout << "Failed to open gnuplot file for mesh plot... \n"; exit(1);}
	gnuplot_file << "set terminal png\n";							// output to "png" file
	gnuplot_file << "set output \"" << polyname << "_mesh.png\"\n";
	gnuplot_file << "set size ratio 0.5\n";
	gnuplot_file << "plot \"" << filename << "\" with lines\n";
	gnuplot_file.close();

	/* Call gnuplot to make plot */
//	string cmd = "gnuplot " + meshfile_name + " --persist";
	string cmd = "gnuplot " + meshfile_name + " \n";
	system(cmd.c_str());
}
