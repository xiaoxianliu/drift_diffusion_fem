#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

#include "mesh.hpp"


using namespace std;

int polyfileSmartInterface(const string polyname, const int num_squares)
{

	/* Nodes */
	double x[4*num_squares + 6], y[4*num_squares+6];		// coordinates of all the nodes
	int num_node_marker=0;
	int* node_marker;

	double x_mid_left = -0.2;					// x-coordinates of the leftmost node on interface
	double x_mid_right = 0.3;					// x-coordinates of the rightmost node on interface

	x[0] = x_mid_left;	y[0] = 0.0;					// vertices define the outer boundary
	x[1] = 1.0;		y[1] = 0.0;
	x[2] = 1.0;		y[2] = 1.0;
	x[3] = x_mid_left;	y[3] = 1.0;
	x[4] = -1.0;		y[4] = 1.0;
	x[5] = -1.0;		y[5] = 0.0;
	
	double delta_y = 1./(2*num_squares + 1);

	for (int i=0; i<num_squares; i++)
	{	x[6+4*i] = x_mid_left;		y[6+4*i] = (2*i+1)*delta_y;	
		x[6+4*i+1] = x_mid_right; 	y[6+4*i+1] = (2*i+1)*delta_y;
		x[6+4*i+2] = x_mid_right;	y[6+4*i+2] = (2*i+2)*delta_y;
		x[6+4*i+3] = x_mid_left;	y[6+4*i+3] = (2*i+2)*delta_y;
	}

	/*Segments and segment attributes */
	int num_seg = 6 + (2*num_squares+1) + 2*num_squares;	// outter segments + interface segments in x direction + interface segments in y direction
	int num_seg_attr = 1;
	int segments[2*num_seg];
	int segment_attr[num_seg];

	for (int i=0; i<5; i++)					// (1) outter segmenets
	{	segments[2*i] = i;	
		segments[2*i+1] = i+1;	
	}
	segments[2*5] = 5;	segments[2*5+1] = 0;

	segment_attr[0] = 2;	segment_attr[1] = 3;	segment_attr[2] = 4;
	segment_attr[3] = 4;	segment_attr[4] = 1;	segment_attr[5] = 2;

	int interface[4*num_squares+2];				// (2) indices of interface nodes
	interface[0] = 0;	interface[4*num_squares+1] = 3;
	for (int i=1;	i<(4*num_squares+1); i++)
	{	interface[i] = 5+i;
	}

	for (int i=0; i<(4*num_squares+1); i++)			//    segments on interface
	{	segments[2*6+2*i] = interface[i];
		segments[2*6+2*i+1] = interface[i+1];
		segment_attr[6+i] = 5;				//    segmenet attributes on interface: all "5"
	}

	/* Regional attributes */
	int num_region_attr = 2;
	double region_x[2]={-0.99, 0.99};
	double region_y[2]={0.01, 0.01};
	int region_attr[2] = {1,2};


	/* Write .poly file */
	WritePolyfile(	polyname,
			4*num_squares+6, x, y, num_node_marker, node_marker,		// vertices
			num_seg, segments, num_seg_attr, segment_attr,			// segments
			num_region_attr, region_x, region_y, region_attr);	
	return 0;
}
