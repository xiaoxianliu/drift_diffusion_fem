#include <cstdlib>
#include <iostream>
#include <vector>
#include "mesh.hpp"
using namespace std;


/*****************************************************************************/
///////////////////////////////////////////////////////////////////////////////////


/* pre-declare each functions */
void ComputeTopology2to0(MeshData &mesh);
void ComputeTopology2to1(MeshData &mesh);
void ComputeTopology2to2(MeshData &mesh);

void ComputeTopology1to0(MeshData &mesh);
void ComputeTopology1to1(MeshData &mesh);
void ComputeTopology1to2(MeshData &mesh);

void ComputeTopology0to0(MeshData &mesh);
void ComputeTopology0to1(MeshData &mesh);
void ComputeTopology0to2(MeshData &mesh);


int ComputeTopology (MeshData &mesh){
	mesh.topology2to0 = 0;
	mesh.topology2to1 = 0;
	mesh.topology2to2 = 0;
	mesh.topology1to0 = 0;
	mesh.topology1to1 = 0;
	mesh.topology1to2 = 0;
	mesh.topology0to0 = 0;
	mesh.topology0to1 = 0;
	mesh.topology0to2 = 0;

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
	// allocate memory for "mesh.topology2to0" array
	mesh.topology2to0 = new vector<int>[mesh.num_elements];

	// Set values to "mesh.topology2to0" array (a copy of "mesh.elements")
	int N = mesh.num_nodes_per_ele;
	for (int i=0; i<mesh.num_elements; i++)
	{	for (int j=0; j<N; j++)
		{	mesh.topology2to0[i].push_back(mesh.elements[i*N+j]);
		}
	}
}

void ComputeTopology1to0 (MeshData &mesh){
	// allocate memory for "mesh.topology1to0" array
	mesh.topology1to0 = new vector<int>[mesh.num_edges];
	// Set values to "mesh.topology1to0" array (a copy of "mesh.edges")
	for (int i=0; i<mesh.num_edges; i++)
	{	for (int j=0; j<2; j++)
		{	mesh.topology1to0[i].push_back(mesh.edges[2*i+j]);		// Only 2 endpoints for each "edge"(i)
		}
	}
}

/* This may include midpoint of edges if elements are not linear */
void ComputeTopology0to2(MeshData &mesh){
	// allocate memory
	mesh.topology0to2 = new vector<int>[mesh.num_nodes];
	// set values

	int N = mesh.num_nodes_per_ele;
	for (int i=0; i<mesh.num_elements; i++)
	{
		for (int j=0; j<N; j++)
		{	int v = mesh.elements[i*N+j];
			mesh.topology0to2[v].push_back(i);
		}
	}
}

void ComputeTopology0to1(MeshData &mesh){
	mesh.topology0to1 = new vector<int>[mesh.num_nodes];				//allocate memory
	for (int i=0; i<mesh.num_edges; i++)						//set values
	{	for (int j=0; j<2; j++)
		{	int v = mesh.edges[2*i+j];
			mesh.topology0to1[v].push_back(i);
		}
	}
}

void ComputeTopology0to0(MeshData &mesh){
	mesh.topology0to0 = new vector<int>[mesh.num_nodes];
	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[2*i];						// identify nodes v0 and v1 of edge "i"
		int v1 = mesh.edges[2*i+1];

		mesh.topology0to0[v0].push_back(v1);
		mesh.topology0to0[v1].push_back(v0);
	}
}


void ComputeTopology1to2(MeshData &mesh){
	if (mesh.topology0to2==0)
		{cout << "Need to compute mesh.topology0to2 before mesh.topology1to2\n";	exit(1);}
	mesh.topology1to2 = new vector<int>[mesh.num_edges];
	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[2*i];
		int v1 = mesh.edges[2*i+1];

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
	if (mesh.topology0to2==0)
		{cout << "Need to compute mesh.topology0to2 before mesh.topology2to1\n";	exit(1);}
	mesh.topology2to1 = new vector<int>[mesh.num_elements];
	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[2*i];
		int v1 = mesh.edges[2*i+1];

		vector<int> T0 = mesh.topology0to2[v0];
		vector<int> T1 = mesh.topology0to2[v1];

		vector<int> T = intersect_vectors(T0, T1);

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
	if (mesh.topology0to2==0)
		{cout << "Need to compute mesh.topology0to2 before mesh.topology2to2\n";	exit(1);}
	mesh.topology2to2 = new vector<int>[mesh.num_elements];
	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[2*i];
		int v1 = mesh.edges[2*i+1];

		vector<int> T0 = mesh.topology0to2[v0];
		vector<int> T1 = mesh.topology0to2[v1];

		vector<int> T = intersect_vectors(T0, T1);

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
	if (mesh.topology0to1==0)
		{cout<< "Need to compute \"topology0to1\" before \"topology1to1\"\n"; exit(1);}
	mesh.topology1to1 = new vector<int> [mesh.num_edges];
	for (int i=0; i<mesh.num_edges; i++)
	{	int v0 = mesh.edges[2*i];
		int v1 = mesh.edges[2*i+1];

		vector<int> E0 = mesh.topology0to1[v0];
		vector<int> E1 = mesh.topology0to1[v1];

		vector<int> E = merge_vectors(E0, E1);
		subtract_one_entry(E, i);

		mesh.topology1to1[i] = E;
	}
}



















