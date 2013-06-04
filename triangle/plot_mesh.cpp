#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include "mesh.hpp"
using namespace std;
using namespace MeshNamespace;


void MeshNamespace::gnuplot_mesh(const MeshData &mesh, const string& filename){

	/* Write data file */
	ofstream output;
	string meshdat_filename = filename + "_mesh.dat";
	output.open(meshdat_filename.c_str());

	if (!output.is_open())
		{cout << "Fail to open "<<filename<<endl; exit(1);}

	for (int i=0; i<mesh.num_elements; i++)
	{	int v0 = mesh.elements[i][0];
		int v1 = mesh.elements[i][1];
		int v2 = mesh.elements[i][2];

		output<<mesh.nodes[v0][0] << "\t" << mesh.nodes[v0][1] << "\n";
		output<<mesh.nodes[v1][0] << "\t" << mesh.nodes[v1][1] << "\n";
		output<<mesh.nodes[v2][0] << "\t" << mesh.nodes[v2][1] << "\n";
		output<<mesh.nodes[v0][0] << "\t" << mesh.nodes[v0][1] << "\n\n";
	}

	output.close();

	/* Write GNUplot command */
	string gnuplot_filename;
	gnuplot_filename = filename + "_mesh.gnuplot";
	ofstream gnuplot_file;
	gnuplot_file.open(gnuplot_filename.c_str());
	if (!gnuplot_file.is_open())								// exit if failed to open file
		{cout << "Failed to open gnuplot file for mesh plot... \n"; exit(1);}
	gnuplot_file << "set terminal png\n";							// output to "png" file
	gnuplot_file << "set output \"" << filename << "_mesh.png\"\n";
	gnuplot_file << "set size ratio 0.5\n";
	gnuplot_file << "plot \"" << meshdat_filename << "\" with lines\n";
	gnuplot_file.close();

	/* Call gnuplot to make plot */
//	string cmd = "gnuplot " + gnuplot_filename + " --persist";
	string cmd = "gnuplot " + gnuplot_filename + " \n";
	system(cmd.c_str());
}




void MeshNamespace::gnuplot_interface(const MeshData &mesh, const string& filename)
{
/*
	vector<int> interface_edge_indices;
	vector<int> interface_node_indices;
	extractInterface(	const MeshData& mesh, \
				vector<int>& interface_edge_indices, \
				vector<int>& interface_node_indices);
*/

	/* 1. Write ".dat" file recording boundary edges and interface edges (indicated by non-zero edge markers) */
	string interface_dat_filename = filename + "_interface.dat";
	ofstream output;
	output.open(interface_dat_filename.c_str());
	if (!output.is_open())
		{std::cout << "Failed to open " << interface_dat_filename << "... Exit!\n";	exit(1);}

	for (int e=0; e<mesh.num_edges; e++)
	{	if (mesh.edge_markers[e]!=0)
		{	int v0 = mesh.topology1to0[e][0];
			int v1 = mesh.topology1to0[e][1];

			output<< mesh.nodes[v0][0] << "\t" << mesh.nodes[v0][1] << "\n";
			output<< mesh.nodes[v1][0] << "\t" << mesh.nodes[v1][1] << "\n\n";
		}
	}

	output.close();


	/* 2. Write gnuplot file */
	string gnuplot_filename = filename + "_interface.gnuplot";
	ofstream gnuplot_file;
	gnuplot_file.open(gnuplot_filename.c_str());
	if (!gnuplot_file.is_open())
		{std::cout << "Failed to open " << gnuplot_filename << ".. Exit!\n";	exit(1);}

	gnuplot_file << "set terminal png\n";
	gnuplot_file << "set output \"" << filename << "_interface.png\"\n";
	gnuplot_file << "set size ratio 0.5\n";
	gnuplot_file << "plot \"" << interface_dat_filename << "\" with lines\n";

	gnuplot_file.close();

	/* 3. Gnuplot */
	string cmd = "gnuplot " + gnuplot_filename + " --persist\n";
	system(cmd.c_str());

}
