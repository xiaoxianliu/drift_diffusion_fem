#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>

using namespace std;

/************************************************************************************/
/*This function write mesh information into a .poly file as input for "Triangle" to generate mesh
/*Note: a few assumptions are made for this function 
/*(1) no node attribute; ("node marker" can exist, though)
/*(2) no hole.
/*(3) no maximum area constraint; it would be stored in a ".area" file if needed
/***********************************************************************************/

int writePolyfile(string filename,
		vector< vector<double> > nodes, vector<int> node_markers,		// vertices
		vector< vector<int> > segments, vector<int> segment_markers,		// segments
		vector< vector<double> > regions, vector<int> region_markers)		// regional attributes
{
	if (nodes.size()==0 || segments.size()==0)
		{cout << "Invalid input: vertex/segmenet number have to be non-zero!\n";	return 1; }

	/* Open .poly file */
	string polyname = filename + ".poly";						// by default, filename has no extension
	ofstream polyfile;
	polyfile.open(polyname.c_str());
	if (!polyfile.is_open())
		{cout << "Failed to open " << polyname << "\n"; exit(1);}

	/* 1. Write vertices */
	int num_marker_per_node = 0;
	if (node_markers.size()>0)
	{	if(node_markers.size() != nodes.size())	
		{cout << "Size of node marker is " << node_markers.size() << "\n";
		 cout << "Size of nodes is "<< nodes.size() << "\n";
		 cout << "They have to be equal!\n";	exit(1);
		}
		else num_marker_per_node=1;
	}
		

	polyfile << "#num of vertices is " << nodes.size() << "; dimension is 2; " \
		  << "num of attributes is 0; " << "num of markers is " << num_marker_per_node << "\n";	// no "attributes" for node

	polyfile << nodes.size() << "\t" << 2 << "\t" << 0 << "\t" << num_marker_per_node << "\n";
	for (int i=0; i< nodes.size(); i++)
	{	polyfile << i << "\t" << nodes[i][0] << "\t" << nodes[i][1];
		if (num_marker_per_node==1)
			polyfile << "\t" << node_markers[i];
		polyfile << "\n";
	}
	polyfile << "\n\n";

	/* 2. Write segments */
	int num_marker_per_segment = 0;
	if ( segment_markers.size()>0 )
		if (segments.size() != segment_markers.size())
		{ cout << "Size of segment markers is " << segment_markers.size() << "\n";
		  cout << "Size of segments is " << segments.size() << "\n";
		  cout << "They have to be equal!\n";	exit(1);
		}
		else num_marker_per_segment = 1;

	polyfile << "# " << segments.size() << " segments; " << num_marker_per_segment << " segment attribute\n";
	polyfile << segments.size() << "\t" << num_marker_per_segment << "\n";
	for (int i=0; i<segments.size(); i++)
	{	polyfile << i << "\t" << segments[i][0] << "\t" << segments[i][1];
		if (num_marker_per_segment == 1)
			polyfile << "\t" << segment_markers[i];
		polyfile << "\n";
	}
	polyfile << "\n\n";




	// Write holes							   // Assume no holes are present
	polyfile << "# 0 hole\n";
	polyfile << 0 << "\n\n\n";




	// Write regional attributes
	int num_marker_per_region = 0;
	if (region_markers.size()>0)
		if (regions.size()!= region_markers.size())
		{ cout << "Size of region markers is "<<region_markers.size() << "\n";
		  cout << "Number of regions is " << regions.size() << "\n";
		  cout << "They have to be equal!\n";	exit(1);
		}
		else	num_marker_per_region = 1;

	polyfile << "# " << regions.size() << " regions\n"; 
	polyfile << regions.size() << "\n";
	for (int i=0; i<regions.size(); i++)
	{	polyfile << i << "\t" << regions[i][0] << "\t" << regions[i][1];
		if ( num_marker_per_region > 0 )
			polyfile << "\t" << region_markers[i] << "\t" << -1 ;	// no maximum area constrait (set to "-1")
		polyfile << "\n";
	}

	polyfile.close();

	return 0;
}
