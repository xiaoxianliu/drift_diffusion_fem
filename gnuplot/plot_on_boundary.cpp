#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include <armadillo>

#include "../triangle/mesh.hpp"

namespace my_gnuplot
{


/* Plot a arma::vec defined over "Omega" only on a segment marked by "seg_marker" */
int plot_ArmaVec_on_DirichletBoundary(	const my_mesh::MeshData &mesh,
					int boundary_marker,
					const arma::vec &u, 
					std::string filename, 
					std::string terminal)
{
using namespace std;
	// 0. Check boundary_marker and identify the corresponding boundary nodes
	if (boundary_marker !=1 && boundary_marker !=3)
	{	std::cout << "Input boundary marker is " << boundary_marker << "; "
			<< "it has to be either 1 or 3 for Dirichlet boundaries.\n";
		exit(1);
	}
	std::vector<int> boundary_nodes;
	if (boundary_marker == 1)	boundary_nodes = mesh.boundary1_nodes;
	else	boundary_nodes = mesh.boundary3_nodes;

	// 1. Write ".dat" file 
	// 1.1 Obtain ".dat" filename
	stringstream dat_filename_stream;
	dat_filename_stream << filename << "(on_boundary_" << boundary_marker << ").dat";
	string dat_filename = dat_filename_stream.str();
	// 1.2 Write to ".dat" file
	ofstream dat_file;
	dat_file.open(dat_filename.c_str());
	if (!dat_file.is_open())
	{	cout << "Failed to open " << dat_filename << "\n";	exit(1);	}

	for (int i=0; i<boundary_nodes.size(); i++)
	{
		int node_index = boundary_nodes[i];
		double y = mesh.nodes[node_index][1];
		dat_file << y << "\t" << u(node_index) << "\n";
	}
	dat_file.close();


	// 2. Write ".gnuplot" file 
	// 2.1 Form ".gnuplot" file stream
	stringstream gnuplot_filename_stream;
	gnuplot_filename_stream << filename << "(on_boundary_" << boundary_marker << ").gnuplot";
	string gnuplot_filename = gnuplot_filename_stream.str();
	// 2.2 Write gnuplot file
	ofstream gnuplot_file;
	gnuplot_file.open(gnuplot_filename.c_str());
	if (!gnuplot_file.is_open())
	{	cout << "Failed to open " << gnuplot_filename << "\n";	exit(1);}
	

	if (terminal == "wxt")
	{	gnuplot_file << "set terminal wxt\n"
				<< "plot \"" << dat_filename << "\" with lines\n";
	}
	else if (terminal == "png")
	{	gnuplot_file << "set terminal png\n"
			 << "set output \"" << filename << "(on_boundary_" << boundary_marker << ").png\"\n"
			 << "plot \"" << dat_filename << "\" with lines\n";
	}
	gnuplot_file.close();

	/* Run gnuplot */
	string cmd = "gnuplot \"" + gnuplot_filename + "\" --persist\n";
	system(cmd.c_str());

	return 0;
}





/************************************************************************************************************************************/





/* Plot arma::vec that's only defined on the some boundary */
int plot_BoundaryArmaVec_on_DirichletBoundary(	const my_mesh::MeshData &mesh,
						int boundary_marker, 
						const arma::vec &u, 
						std::string filename, 
						std::string terminal)
{
// Note: entries of "arma::vec u" MUST BE consistent with the order of "mesh.interface_nodes" 

using namespace std;
	// 0. Check boundary marker and identify correponding boundary nodes
	if (boundary_marker !=1 && boundary_marker !=3)
	{	std::cout << "Input boundary marker is " << boundary_marker << "; "
			<< "it has to be either 1 or 3 for Dirichlet boundaries.\n";
		exit(1);
	}
	std::vector<int> boundary_nodes;
	if (boundary_marker == 1)	boundary_nodes = mesh.boundary1_nodes;
	else	boundary_nodes = mesh.boundary3_nodes;


	// 1. Write ".dat" file */
	// 1.1 Obtain ".dat" filename
	stringstream dat_filename_stream;
	dat_filename_stream << filename << "(on_boundary_" << boundary_marker << ").dat";
	string dat_filename = dat_filename_stream.str();
	// 1.2 Write to ".dat" file
	ofstream dat_file;
	dat_file.open(dat_filename.c_str());
	if (!dat_file.is_open())
	{	cout << "Failed to open " << dat_filename << "\n";	exit(1);	}

	for (int i=0; i<boundary_nodes.size(); i++)
	{
		int node_index = boundary_nodes[i];
		double y = mesh.nodes[node_index][1];
		dat_file << y << "\t" << u(i) << "\n";
	}
	dat_file.close();


	// 2. Write ".gnuplot" file */
	// 2.1 Form ".gnuplot" file stream
	stringstream gnuplot_filename_stream;
	gnuplot_filename_stream << filename << "(on_boundary_" << boundary_marker << ").gnuplot";
	string gnuplot_filename = gnuplot_filename_stream.str();
	// 2.2 Write gnuplot file
	ofstream gnuplot_file;
	gnuplot_file.open(gnuplot_filename.c_str());
	if (!gnuplot_file.is_open())
	{	cout << "Failed to open " << gnuplot_filename << "\n";	exit(1);}


	if (terminal == "wxt")
	{	gnuplot_file << "set terminal wxt\n";
		gnuplot_file << "plot \"" << dat_filename << "\" with lines\n";
	}
	else if (terminal == "png")
	{	gnuplot_file << "set terminal png\n"
			 << "set output \"" << filename << "(on_boundary_" << boundary_marker << ").png\"\n"
			 << "plot \"" << dat_filename << "\" with lines\n";
	}
	gnuplot_file.close();

	/* Run gnuplot */
	string cmd = "gnuplot \"" + gnuplot_filename + "\" --persist\n";
	system(cmd.c_str());

	return 0;
}




















}
