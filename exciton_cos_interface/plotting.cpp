#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

#include <armadillo>

#include "../triangle/mesh.hpp"


/* Plot arma::vec over the entire domain "Omega" */
int plot_ArmaVec(const my_mesh::MeshData& mesh, const arma::vec u, const std::string& filename)
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

	gnuplot_fstream	<< "set terminal png\n"
			<< "set output \"" << filename << "_solution.png\"\n"
			<< "splot \"" << output_name << "\" with lines\n";

	gnuplot_fstream.close();


	/* Run GNUplot file */
	string cmd = "gnuplot " + gnuplot_filename + " --persist\n";
	system(cmd.c_str());

	return 0;
}






/* Plot a arma::vec defined over "Omega" only on the interface "Gamma" */
int plot_ArmaVec_on_Interface(const my_mesh::MeshData &mesh, const arma::vec &u, std::string filename)
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
	
	gnuplot_file << "set terminal png\n";
	gnuplot_file << "set output \"" << filename << "(on_interface).png\"\n";
	gnuplot_file << "plot \"" << dat_filename << "\" with lines\n";

	gnuplot_file.close();

	/* Run gnuplot */
	string cmd = "gnuplot \"" + gnuplot_filename + "\" --persist\n";
	system(cmd.c_str());

	return 0;
}











/* Plot arma::vec that's only defined on the interface Gamma */
int plot_InterfaceArmaVec(const my_mesh::MeshData &mesh, const arma::vec &u, std::string filename)
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
	
	gnuplot_file << "set terminal png\n";
	gnuplot_file << "set output \"" << filename << "(on_interface).png\"\n";
	gnuplot_file << "plot \"" << dat_filename << "\" with lines\n";

	gnuplot_file.close();

	/* Run gnuplot */
	string cmd = "gnuplot \"" + gnuplot_filename + "\" --persist\n";
	system(cmd.c_str());

	return 0;
}









/* Plot std::vector that's only defined on the interface Gamma */
int plot_InterfaceSTLVector(const my_mesh::MeshData &mesh, const std::vector<double> &u, std::string filename)
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
	
	gnuplot_file << "set terminal png\n";
	gnuplot_file << "set output \"" << filename << "(on_interface).png\"\n";
	gnuplot_file << "plot \"" << dat_filename << "\" with lines\n";

	gnuplot_file.close();

	/* Run gnuplot */
	string cmd = "gnuplot \"" + gnuplot_filename + "\" --persist\n";
	system(cmd.c_str());

	return 0;
}
