#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

#include <armadillo>

#include "../triangle/mesh.hpp"

namespace my_gnuplot
{


/* Plot a arma::vec defined over "Omega" only on the interface "Gamma" */
int plot_ArmaVec_on_Interface(const my_mesh::MeshData &mesh, const arma::vec &u, std::string filename, std::string terminal)
{
using namespace std;

	/* Write ".dat" file */
	string dat_filename = filename + "(on_interface).dat";
	ofstream dat_file;
	dat_file.open(dat_filename.c_str());
	if (!dat_file.is_open())
	{	cout << "Failed to open " << dat_filename << " in function \"plotSolution_on_Interface\"\n";	exit(1);	}

	for (int i=0; i<mesh.interface_nodes.size(); i++)
	{
		int node_index = mesh.interface_nodes[i];
		double y = mesh.nodes[node_index][1];
		dat_file << y << "\t" << u(node_index) << "\n";
	}
	dat_file.close();


	/* Write ".gnuplot" file */
	string gnuplot_filename = filename + "(on_interface).gnuplot";
	ofstream gnuplot_file;
	gnuplot_file.open(gnuplot_filename.c_str());
	if (!gnuplot_file.is_open())
	{	cout << "Failed to open " << gnuplot_filename << " in function \"plotSolution_on_Interface\"\n";	exit(1);}
	

	if (terminal == "wxt")
	{	gnuplot_file << "set terminal wxt\n"
				<< "plot \"" << dat_filename << "\" with lines\n";
	}
	else if (terminal == "png")
	{	gnuplot_file << "set terminal png\n";
		gnuplot_file << "set output \"" << filename << "(on_interface).png\"\n";
		gnuplot_file << "plot \"" << dat_filename << "\" with lines\n";
	}
	gnuplot_file.close();

	/* Run gnuplot */
	string cmd = "gnuplot \"" + gnuplot_filename + "\" --persist\n";
	system(cmd.c_str());

	return 0;
}











/* Plot arma::vec that's only defined on the interface Gamma */
int plot_InterfaceArmaVec(const my_mesh::MeshData &mesh, const arma::vec &u, std::string filename, std::string terminal)
{
/* Note: "arma::vec u" MUST BE consistent with the order of "mesh.interface_nodes" */

using namespace std;

	/* Write ".dat" file */
	string dat_filename = filename + "(on_interface).dat";
	ofstream dat_file;
	dat_file.open(dat_filename.c_str());
	if (!dat_file.is_open())
	{	cout << "Failed to open " << dat_filename << " in function \"plot_InterfaceArmaVec\"\n";	exit(1);	}

	for (int i=0; i<mesh.interface_nodes.size(); i++)
	{
		int node_index = mesh.interface_nodes[i];
		double y = mesh.nodes[node_index][1];
		dat_file << y << "\t" << u(i) << "\n";
	}
	dat_file.close();


	/* Write ".gnuplot" file */
	string gnuplot_filename = filename + "(on_interface).gnuplot";
	ofstream gnuplot_file;
	gnuplot_file.open(gnuplot_filename.c_str());
	if (!gnuplot_file.is_open())
	{	cout << "Failed to open " << gnuplot_filename << " in function \"plot_InterfaceArmaVec\"\n";	exit(1);}


	if (terminal == "wxt")
	{	gnuplot_file << "set terminal wxt\n";
		gnuplot_file << "plot \"" << dat_filename << "\" with lines\n";
	}
	else if (terminal == "png")
	{	gnuplot_file << "set terminal png\n";
		gnuplot_file << "set output \"" << filename << "(on_interface).png\"\n";
		gnuplot_file << "plot \"" << dat_filename << "\" with lines\n";
	}
	gnuplot_file.close();

	/* Run gnuplot */
	string cmd = "gnuplot \"" + gnuplot_filename + "\" --persist\n";
	system(cmd.c_str());

	return 0;
}









/* Plot std::vector that's only defined on the interface Gamma */
int plot_InterfaceSTLVector(const my_mesh::MeshData &mesh, const std::vector<double> &u, std::string filename, std::string terminal)
{
/* Note: "std::vector<double> u" MUST BE consistent with the order of "mesh.interface_nodes" */

using namespace std;

	/* Write ".dat" file */
	string dat_filename = filename + "(on_interface).dat";
	ofstream dat_file;
	dat_file.open(dat_filename.c_str());
	if (!dat_file.is_open())
	{	cout << "Failed to open " << dat_filename << " in function \"plotSolution_on_Interface\"\n";	exit(1);	}

	for (int i=0; i<mesh.interface_nodes.size(); i++)
	{
		int node_index = mesh.interface_nodes[i];
		double y = mesh.nodes[node_index][1];
		dat_file << y << "\t" << u[i] << "\n";
	}
	dat_file.close();


	/* Write ".gnuplot" file */
	string gnuplot_filename = filename + "(on_interface).gnuplot";
	ofstream gnuplot_file;
	gnuplot_file.open(gnuplot_filename.c_str());
	if (!gnuplot_file.is_open())
	{	cout << "Failed to open " << gnuplot_filename << " in function \"plotSolution_on_Interface\"\n";	exit(1);}


	if (terminal == "wxt")
	{	gnuplot_file << "set terminal wxt\n";
		gnuplot_file << "plot \"" << dat_filename << "\" with lines\n";
	}
	else if (terminal == "png")
	{	gnuplot_file << "set terminal png\n";
		gnuplot_file << "set output \"" << filename << "(on_interface).png\"\n";
		gnuplot_file << "plot \"" << dat_filename << "\" with lines\n";
	}

	gnuplot_file.close();

	/* Run gnuplot */
	string cmd = "gnuplot \"" + gnuplot_filename + "\" --persist\n";
	system(cmd.c_str());

	return 0;
}


















}
