#include <iostream>
#include <string>

#include <armadillo>
#include "mesh.hpp"		// In "/home/xiaoxian/programs_cpp/triangle/

#include "exciton_diffusion.hpp"

using namespace std;

int main(int argc, char *argv[])
{
	/*********** Generate mesh */
	/* 1. Read in name of input file for "triangle" */
	string polyname = "rectangle_smart";
	int num_squares_on_interface = 0;
	double drainage_on_interface;

	Solve(polyname, num_squares_on_interface, drainage_on_interface);

	return 0;
}
