#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>

#include "mesh.hpp"

int ReadElements(MeshData &mesh, string ele_file_name){
	ifstream input_file;
	input_file.open(ele_file_name.c_str());
	if (!input_file.is_open())
		{cout << "Fail to open .ele file " << ele_file_name << endl;	exit(1);	}

	/* 1. Read first line of basic information of elements */
	string line;
	getline(input_file, line);
	stringstream linestream(line);
	int num_elements, num_nodes_per_ele, num_attributes_per_ele;
	linestream >> num_elements >> num_nodes_per_ele >> num_attributes_per_ele;
	if ( (num_attributes_per_ele!=0) and (num_attributes_per_ele!=1) )
		{cout << "# of element attributes has to be either 0 or 1 (assumbed to be used as marker) " << endl;	exit(1);}

	mesh.num_elements = num_elements;
	mesh.num_nodes_per_ele = num_nodes_per_ele;
	mesh.num_attributes_per_ele = num_attributes_per_ele;


	/* 2. Read in element information */
	while (input_file.good())
	{	getline(input_file, line);

		/* Skip irrelevant lines */
		if (line.length()==0)
			{continue;}
		if (line[0]=='#')
			{continue;}

		/* Read information to temporary variables */
		stringstream linestream(line);
		int index;
		vector<int> ele(num_nodes_per_ele);
		int marker=0;

		linestream >> index;							// read index
		for (int i=0; i<num_nodes_per_ele; i++)					// read nodal index
			{ linestream >> ele[i];	}
		if (num_attributes_per_ele == 1)					// read element marker
			{linestream >> marker;	}


		/* Finally, set relevant elements in "mesh" to the value of the temporary variables above */
		mesh.elements.push_back(ele);
		mesh.element_markers.push_back(marker);

/* Debug
		cout << "index " << index << ": ";
		for (int i=0; i<num_nodes_per_ele; i++)
			{ cout << mesh.elements[index][i] << ", ";	}
		cout << mesh.element_markers[index] << endl;
*/
	}

	input_file.close();

	/* Final check on total number of elements */
	if (mesh.num_elements != mesh.elements.size())
		{  cout << "\"mesh.num_elements\" is " << mesh.num_elements << "\n";
		   cout << "\"mesh.elements\" is a cpp vector of size " << mesh.elements.size() << "\n";
		   cout << "They have to be equal! Exit... \n";
		   exit(1);
		}
	return 0;
}
