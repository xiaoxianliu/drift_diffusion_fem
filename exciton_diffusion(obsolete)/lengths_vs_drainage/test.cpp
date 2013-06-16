#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <armadillo>
#include "mesh.hpp"		// In "/home/xiaoxian/programs_cpp/triangle/

#include "exciton_diffusion.hpp"

using namespace std;

int main(int argc, char *argv[])
{
	int max_num_squares=1;

	double interface_lengths[max_num_squares];
	double drainage[max_num_squares];

	for (int i=0; i<(max_num_squares+1); i++)
	{	interface_lengths[i] = i + 1.0;
		double drainage_i;

		std::string filename="rectangle_smart";
		std::ostringstream filename_stream;
		filename_stream << filename << i;
		filename = filename_stream.str();

		Solve(filename, i, drainage_i);
		drainage[i] = drainage_i;

	}

	/* Save result in .dat file */
	string output_name = "output";
	string dat_file_name= output_name + ".dat";
	ofstream dat_file_stream;
	dat_file_stream.open(dat_file_name.c_str());

	dat_file_stream << "#interface lengths vs drainage\n";
	for (int i=0; i<(max_num_squares+1); i++)
	{	dat_file_stream << interface_lengths[i] << "\t" << drainage[i] << "\n";
	}

	dat_file_stream.close();

	/* Plot result */
	string gnuplot_file_name= output_name + ".gnuplot";
	ofstream gnuplot_file_stream;
	gnuplot_file_stream.open(gnuplot_file_name.c_str());

	gnuplot_file_stream << "set terminal png\n";
	gnuplot_file_stream << "set output " << ("\"" + output_name + ".png\"") << "\n";
	gnuplot_file_stream << "plot " << ( "\"" + dat_file_name + "\"" ) << " with linespoints\n";

	gnuplot_file_stream.close();

	string cmd = "gnuplot " + gnuplot_file_name + " --persist\n";
	system(cmd.c_str());

	return 0;
}
