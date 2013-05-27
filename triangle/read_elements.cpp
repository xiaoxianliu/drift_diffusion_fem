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

	/* Read first line of basic information of elements */
	string line;
	getline(input_file, line);
	stringstream linestream(line);
	int num_elements, num_nodes_per_ele, num_attributes;
	linestream >> num_elements >> num_nodes_per_ele >> num_attributes;
	if ( (num_attributes!=0) and (num_attributes!=1) )
		{cout << "# of element attributes has to be either 0 or 1 (used as marker) " << endl;	exit(1);}

	/* Allocate memory for relevant elements of "mesh" */
	mesh.num_elements = num_elements;
	mesh.num_nodes_per_ele = num_nodes_per_ele;
	mesh.num_ele_attributes = num_attributes;
	mesh.elements = new int[num_elements * num_nodes_per_ele];
	mesh.element_markers = new int[num_elements];



/*	cout << "number of elements is " << mesh.num_elements << endl;
	cout << "number of nodes per element is " << mesh.num_nodes_per_ele << endl;
	cout << "number of element attributes is " << mesh.num_ele_attributes << endl;
*/

	/* Read in element information */
	while (input_file.good())
	{	getline(input_file, line);

		/* Skip unrelevant lines */
		if (line.length()==0)
			{continue;}
		if (line[0]=='#')
			{continue;}

		/* Read information to temporary variables */
		stringstream linestream(line);
		int index;
		int nodes[num_nodes_per_ele];
		int marker=0;
		linestream >> index;							// read index
		for (int i=0; i<num_nodes_per_ele; i++)					// read nodal index
			{ linestream >> nodes[i];	}
		if (num_attributes>0)							// read element marker
			{linestream >> marker;	}


		/* Finally, set relevant elements in "mesh" to the value of the temporary variables above */
		for (int i=0; i<num_nodes_per_ele; i++)
			{ mesh.elements[index * num_nodes_per_ele + i] = nodes[i];	}
		mesh.element_markers[index] = marker;

/* Debug
		cout << "index " << index << ": ";
		for (int i=0; i<num_nodes_per_ele; i++)
			{ cout << mesh.elements[index*num_nodes_per_ele + i] << ", ";	}
		cout << mesh.element_markers[index] << endl;
*/


	}

	input_file.close();
	return 0;
}
