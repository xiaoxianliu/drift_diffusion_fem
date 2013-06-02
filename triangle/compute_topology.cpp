#include <cstdlib>
#include <iostream>
#include <vector>
#include "mesh.hpp"


using namespace std;			// std::vector, std::cout, std::endl
using namespace MeshNamespace;

/*****************************************************************************/



/* pre-declare each local functions */
void ComputeTopology2to0(MeshData &mesh);
void ComputeTopology2to1(MeshData &mesh);
void ComputeTopology2to2(MeshData &mesh);

void ComputeTopology1to0(MeshData &mesh);
void ComputeTopology1to1(MeshData &mesh);
void ComputeTopology1to2(MeshData &mesh);

void ComputeTopology0to0(MeshData &mesh);
void ComputeTopology0to1(MeshData &mesh);
void ComputeTopology0to2(MeshData &mesh);


/* Main function to compute every topology */
int MeshNamespace::ComputeTopology (MeshData &mesh){

	ComputeTopology2to0(mesh);
	ComputeTopology1to0(mesh);
	ComputeTopology0to2(mesh);
	ComputeTopology0to1(mesh);
	ComputeTopology0to0(mesh);

	/* The 3 topologies below need information of "mesh.topology0to2" before being computed */
	ComputeTopology1to2(mesh);
	ComputeTopology2to1(mesh);
	ComputeTopology2to2(mesh);

	/* Needs information of "mesh.topology1to1" */
	ComputeTopology1to1(mesh);

	return 0;
}










/***********************************************************************************/
/*************** Define functions computing each topology **************************/


/* This may include midpoint of edges if elements are not linear */
void ComputeTopology2to0 (MeshData &mesh){
	mesh.topology2to0.clear();		// clear potential existing data

	// A copy of "mesh.elements"
	mesh.topology2to0.assign(mesh.elements.begin(), mesh.elements.end());
}

void ComputeTopology1to0 (MeshData &mesh){
	mesh.topology1to0.clear();		// clear potential existing data

	// Set values to "mesh.topology1to0" array (a copy of "mesh.edges")
	mesh.topology1to0.assign(mesh.edges.begin(), mesh.edges.end());

}

/* This may include midpoint of edges if elements are not linear */
void ComputeTopology0to2(MeshData &mesh){

	mesh.topology0to2.clear();
	mesh.topology0to2.resize(mesh.num_nodes);

	for (int i=0; i<mesh.elements.size(); i++)
	{
		for (int j=0; j<mesh.elements[i].size(); j++)
		{	int v = mesh.elements[i][j];
			mesh.topology0to2[v].push_back(i);
		}
	}
}

void ComputeTopology0to1(MeshData &mesh){
	mesh.topology0to1.clear();							//clear and resize "topology0to1"
	mesh.topology0to1.resize(mesh.num_nodes);
	for (int i=0; i<mesh.num_edges; i++)						//set values
	{	for (int j=0; j<2; j++)
		{	int v = mesh.edges[i][j];
			mesh.topology0to1[v].push_back(i);
		}
	}
}

void ComputeTopology0to0(MeshData &mesh){
	mesh.topology0to0.clear();							//clear and resize "topology0to0"
	mesh.topology0to0.resize(mesh.num_nodes);

	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[i][0];						// identify nodes v0 and v1 of edge "i"
		int v1 = mesh.edges[i][1];

		mesh.topology0to0[v0].push_back(v1);
		mesh.topology0to0[v1].push_back(v0);
	}
}


void ComputeTopology1to2(MeshData &mesh){
	if (mesh.topology0to2.size()==0)
		{cout << "Need to compute mesh.topology0to2 before mesh.topology1to2\n";	exit(1);}

	mesh.topology1to2.clear();							//clear and resize "topology1to2"
	mesh.topology1to2.resize(mesh.num_edges);

	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[i][0];
		int v1 = mesh.edges[i][1];

		vector<int> T0 = mesh.topology0to2[v0];
		vector<int> T1 = mesh.topology0to2[v1];

		vector<int> T = intersect_vectors(T0, T1);
		if (T.size()!=2 && T.size()!=1)
			{	cout<<"Computing \"topology1to2\"..." << endl;

				cout<<"Nodes of edge "<<i<<" are " << v0 <<","<<v1<<endl;

				cout<<"Neighboring elements of "<<v0<< " are ";
				for (int ii=0; ii<T0.size(); ii++)
					cout << T0[ii]<<" ";
				cout<<endl;

				cout<<"Neighboring elements of "<<v1<< " are ";
				for (int ii=0; ii<T1.size(); ii++)
					cout << T1[ii]<<" ";
				cout<<endl;

				cout<<"The elements shared by "<<v0<<" and "<<v1<<" are:";
				for (int ii=0; ii<T.size(); ii++)
					cout << T[ii]<<" ";
				cout<<endl;

				cout<<"T0 and T1 must share only 1 or 2 numbers... exit!"<<endl;
				exit(1);
			}
		else
			mesh.topology1to2[i]=T;
	}
}

void ComputeTopology2to1(MeshData &mesh){
	if (mesh.topology0to2.size()==0)
		{cout << "Need to compute mesh.topology0to2 before mesh.topology2to1\n";	exit(1);}

	mesh.topology2to1.clear();							// clear and resize "topology2to1"
	mesh.topology2to1.resize(mesh.num_elements);

	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[i][0];
		int v1 = mesh.edges[i][1];

		vector<int> T0 = mesh.topology0to2[v0];
		vector<int> T1 = mesh.topology0to2[v1];

		vector<int> T = intersect_vectors(T0, T1);				// elements shared by v0 and v1

		if (T.size()!=2 && T.size()!=1)
			{	cout<<"Computing \"topology2to1\"... " << endl;

				cout<<"Nodes of edge "<<i<<" are " << v0 <<","<<v1<<endl;

				cout<<"Neighboring elements of "<<v0<< " are ";
				for (int ii=0; ii<T0.size(); ii++)
					cout << T0[ii]<<" ";
				cout<<endl;

				cout<<"Neighboring elements of "<<v1<< " are ";
				for (int ii=0; ii<T1.size(); ii++)
					cout << T1[ii]<<" ";
				cout<<endl;

				cout<<"The elements shared by "<<v0<<" and "<<v1<<" are:";
				for (int ii=0; ii<T.size(); ii++)
					cout << T[ii]<<" ";
				cout<<endl;

				cout<<"T0 and T1 must share only 1 or 2 numbers... exit!"<<endl;
				exit(1);
			}
		else
		{	for (int j=0; j<T.size(); j++)
				mesh.topology2to1[ T[j] ].push_back(i);
		}
	}
}


void ComputeTopology2to2(MeshData &mesh){
	if (mesh.topology0to2.size()==0)
		{cout << "Need to compute mesh.topology0to2 before mesh.topology2to2\n";	exit(1);}

	mesh.topology2to2.clear();							// clear and reszie "topology2to2"
	mesh.topology2to2.resize(mesh.num_elements);

	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[i][0];
		int v1 = mesh.edges[i][1];

		vector<int> T0 = mesh.topology0to2[v0];
		vector<int> T1 = mesh.topology0to2[v1];

		vector<int> T = intersect_vectors(T0, T1);				// elements shared by "v0" and "v1"

		if (T.size()!=2 && T.size()!=1)
			{	cout<<"Computing \"topology2to2\".." << endl;

				cout<<"Nodes of edge "<<i<<" are " << v0 <<","<<v1<<endl;

				cout<<"Neighboring elements of "<<v0<< " are ";
				for (int ii=0; ii<T0.size(); ii++)
					cout << T0[ii]<<" ";
				cout<<endl;

				cout<<"Neighboring elements of "<<v1<< " are ";
				for (int ii=0; ii<T1.size(); ii++)
					cout << T1[ii]<<" ";
				cout<<endl;

				cout<<"The elements shared by "<<v0<<" and "<<v1<<" are:";
				for (int ii=0; ii<T.size(); ii++)
					cout << T[ii]<<" ";
				cout<<endl;

				cout<<"T0 and T1 must share only 1 or 2 numbers... exit!"<<endl;
				exit(1);
			}
		else if (T.size()==2)
			{	int t0 = T[0];
				int t1 = T[1];

				mesh.topology2to2[t0].push_back(t1);
				mesh.topology2to2[t1].push_back(t0);
			}
	}
}



void ComputeTopology1to1(MeshData &mesh){
	if (mesh.topology0to1.size()==0)
		{cout<< "Need to compute \"topology0to1\" before \"topology1to1\"\n"; exit(1);}

	mesh.topology1to1.clear();							// clear and resize "topology1to1"
	mesh.topology1to1.resize(mesh.num_edges);

	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[i][0];
		int v1 = mesh.edges[i][1];

		vector<int> E0 = mesh.topology0to1[v0];
		vector<int> E1 = mesh.topology0to1[v1];

		vector<int> E = merge_vectors(E0, E1);
		subtract_one_entry(E, i);

		mesh.topology1to1[i] = E;
	}
}



















