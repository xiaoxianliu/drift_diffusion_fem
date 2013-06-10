#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include "mesh.hpp"


MeshNamespace::MeshData MeshNamespace::refineMesh(	const MeshNamespace::MeshData &old_mesh,
							std::string filename, 
							int level)
// example: for mesh corresponding to "mesh.1.ele",  "filename" should be "mesh" and "level" should be "1"
{using namespace std;

	/* Form the "filename" which includes current "level" of refinement. */
	stringstream old_mesh_filename_stream;
	old_mesh_filename_stream << filename << "." << level;
	string old_mesh_filename;
	old_mesh_filename = old_mesh_filename_stream.str();
	

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
	string old_ele_filename = old_mesh_filename + ".ele";			// file name of current level of ".ele" file
	ifstream old_ele_file;
	old_ele_file.open(old_ele_filename.c_str());
	if (!old_ele_file.is_open())
	{	cout << "Failed to open " << old_ele_filename << "...Exit!\n"; 	exit(1);	}

	double num_elements_file;
	old_ele_file >> num_elements_file;
	old_ele_file.close();

	// Check if the two "number of elements" are consistent. If not, exit.
	if ( num_elements != num_elements_file )
	{	cout << "Number of elements in MeshData old_mesh is " << num_elements << "\n";
		cout << "Number of elements in \"" << old_ele_filename << " is " << num_elements_file << "\n";
		cout << "The should be equal! \n";
		exit(1);
	}




	/* 1. Write ".area" file according to the refinement info in "old_mesh.refinement_markers" */
	string old_area_filename = old_mesh_filename + ".area";
	ofstream old_area_file;
	old_area_file.open(old_area_filename.c_str());
	if ( !old_area_file.is_open() )
	{	cout << "Failed to open " << old_area_filename << ". Exit...\n";	exit(1);	}

	// Write number of elements to .area file
	old_area_file << num_elements << "\n";

	// Write area constraint on each element if needed
	for (int i=0; i < old_mesh.num_elements; i++)
	{
		int marker = old_mesh.refinement_markers[i];
		switch (marker)
		{
			case 0:
				old_area_file << i << "\t" << -1 << "\n";			// do nothing
				break;
			case 1:
				old_area_file << i << "\t" << 0.8 * old_mesh.ele_areas[i] << "\n";		// smaller than current area
				break;
			default:
				cout << i << "-th marker for refinement is " << marker << "; it should be either 0 or 1!\n";
				exit(1);
		}
	}

	old_area_file.close();


	/* 2. Refine mesh */
	string cmd = "/home/xiaoxian/bin/triangle/triangle -rqepzAa " + old_mesh_filename + "\n";
	system(cmd.c_str());


	/* 3. Process new mesh
	/* 3.1 read in new mesh; 3.2 compute topological and geometric properties; 3.3 plot mesh and interface 
	*/
	
	MeshData new_mesh;

	// Form "new_mesh_filename", e.g. "example.2" when "filename"="example" and "level" = 1.
	stringstream new_mesh_filename_stream;
	new_mesh_filename_stream << filename << "." << level+1 ;
	string new_mesh_filename = new_mesh_filename_stream.str();

	string new_node_file_name= new_mesh_filename + ".node";
	string new_edge_file_name= new_mesh_filename + ".edge";
	string new_ele_file_name = new_mesh_filename + ".ele";

	ReadNodes(new_mesh, new_node_file_name);
	ReadEdges(new_mesh, new_edge_file_name);
	ReadElements(new_mesh, new_ele_file_name);

	ComputeTopology(new_mesh);

	ComputeEdgeLengths(new_mesh);
	ComputeElementAreas(new_mesh);

	gnuplot_mesh(new_mesh, new_mesh_filename);
	gnuplot_interface(new_mesh, new_mesh_filename);



	return new_mesh;
}
