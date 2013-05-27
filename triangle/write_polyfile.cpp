#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

/************************************************************************************/
/*This function write mesh information into a .poly file as input for "Triangle" to generate mesh
/*Note: a few assumptions are made for this function 
/*(1) no node attribute; ("node marker" can exist, though)
/*(2) no hole.
/***********************************************************************************/

int WritePolyfile(string polyname,
		int num_nodes, double x[], double y[], int num_node_marker, int node_marker[],		// vertices
		int num_seg, int segments[], int num_seg_attr, int seg_attr[],				// segments
		int num_region_attr, double region_x[], double region_y[], int region_attr[])		// regional attributes
{
	if (num_nodes==0 || num_seg==0)
		{cout << "Invalid input: vertex/segmenet number have to be non-zero!\n";	return 1; }

	polyname += ".poly";								// by default, polyname has no extension
	ofstream polyfile;
	polyfile.open(polyname.c_str());
	if (!polyfile.is_open())
		{cout << "Failed to open " << polyname << "\n"; return 1;}

	// Write vertices
	polyfile << "#num of vertices is " << num_nodes << "; dimension is 2; " \
		  << "num of attributes is 0; " << "num of markers is " << num_node_marker << "\n";

	polyfile << num_nodes << "\t" << 2 << "\t" << 0 << "\t" << num_node_marker << "\n";
	for (int i=0; i< num_nodes; i++)
	{	polyfile << i << "\t" << x[i] << "\t" << y[i];
		if (num_node_marker==1)
			polyfile << "\t" << node_marker[i];
		polyfile << "\n";
	}
	polyfile << "\n\n";

	// Write segments
	polyfile << "# " << num_seg << " segments; " << num_seg_attr << " segment attribute\n";
	polyfile << num_seg << "\t" << num_seg_attr << "\n";
	for (int i=0; i<num_seg; i++)
	{	polyfile << i << "\t" << segments[2*i] << "\t" << segments[2*i+1];
		if (num_seg_attr == 1)
			polyfile << "\t" << seg_attr[i];
		polyfile << "\n";
	}
	polyfile << "\n\n";

	// Write holes							   // Assume no holes are present
	polyfile << "# 0 hole\n";
	polyfile << 0 << "\n\n\n";


	// Write regional attributes
	polyfile << "# " << num_region_attr << " regional attributes\n"; 
	for (int i=0; i<num_region_attr; i++)
	{	polyfile << i << "\t" << region_x[i] << "\t" << region_y[i] << "\t"
			 << region_attr[i] << "\t" << -1 << "\n";
	}

	polyfile.close();

	return 0;
}
