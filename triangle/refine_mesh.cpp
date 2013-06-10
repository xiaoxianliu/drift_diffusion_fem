#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include "mesh.hpp"


MeshNamespace::MeshData refineMesh (	const MeshNamespace::MeshData &old_mesh,
					std::string filename)			// example: for mesh corresponding to "mesh.1.ele", 
										//  "filename" should be "mesh.1"
{using namespace std;

	/* 0. Pre-checking consistency */
	int num_elements = old_mesh.num_elements;

	/* 0.1 Check if "old_mesh.refinement_markers" is properly set */
	if ( old_mesh.refinement_markers.size() != num_elements )
	{	cout << "Number of refinement markers is " << old_mesh.refinement_markers.size() 
			<< ", whereas number of elements is " << num_elements << ".\n";
		cout << "They should be equal! Exit...\n";
		exit(1);
	}

	/* 0.2 Check if information in "old_mesh" is consistent with information in ".ele" file */
	// Read in "number of elements" from the ".ele" file of given "filename"
	string ele_filename = filename + ".ele";			// file name of current level of ".ele" file
	ifstream ele_file;
	ele_file.open(ele_filename.c_str());
	if (!ele_file.is_open())
	{	cout << "Failed to open " << ele_filename << "...Exit!\n"; 	exit(1);	}

	double num_elements_file;
	ele_file >> num_elements_file;
	ele_file.close();

	// Check if the two "number of elements" are consistent. If not, exit.
	if ( num_elements != num_elements_file )
	{	cout << "Number of elements in MeshData old_mesh is " << num_elements << "\n";
		cout << "Number of elements in \"" << ele_filename << " is " << num_elements_file << "\n";
		cout << "The should be equal! \n";
		exit(1);
	}




	/* 1. Write ".area" file according to the refinement info in "old_mesh.refinement_markers" */
	string area_filename = filename + ".area";
	ofstream area_file;
	area_file.open(area_filename.c_str());
	if ( !area_file.is_open() )
	{	cout << "Failed to open " << area_filename << ". Exit...\n";	exit(1);	}

	// Write number of elements to .area file
	area_file << num_elements << "\n";

	// Write area constraint on each element if needed
	for (int i=0; i < old_mesh.num_elements; i++)
	{
		int marker = old_mesh.refinement_markers[i];
		switch (marker)
		{
			case 0:
				area_file << i << "\t" << -1 << "\n";			// do nothing
				break;
			case 1:
				area_file << i << "\t" << 0.8 * old_mesh.ele_areas[i] << "\n";		// smaller than current area
				break;
			default:
				cout << i << "-th marker for refinement is " << marker << "; it should be either 0 or 1!\n";
				exit(1);
		}
	}

	area_file.close();


	/* 2. Refine mesh */
	string cmd = "/home/xiaoxian/bin/triangle/triangle -qepzAa " + filename + "\n";
	system(cmd.c_str());

}
