#include <vector>
#include "mesh.hpp"

using std::vector;

/* Add an entry(int typed) "n" to given vector<int> v, if n is not in "v" already; otherwise do nothing */
int MeshNamespace::add_one_entry(vector<int> &v, const int &n)
{	bool not_existing = true;

	for (int i=0; i<v.size(); i++)
	{	if (v[i]==n)
			{ not_existing=false;	break;	}
	}
	if (not_existing)
		v.push_back(n);
	return 0;
}

/* Subtract an entry(int typed) "n" from given vector<int> v, if n is in "v"; otherwise, do nothing. */
int MeshNamespace::subtract_one_entry(vector<int> &v, const int &n)
{	for (int i=0; i<v.size(); i++)
		if (v[i]==n)
			v.erase(v.begin()+i);
	return 0;
}

/* Take intersection of two vector<int> object, "v1" and "v2". And return a vector<int> object, "v_out" */
vector<int> MeshNamespace::intersect_vectors(const vector<int> &v1, const vector<int> &v2)
{	vector<int> v_out(v1);
	bool in_both = false;
	for (int i=0; i<v1.size(); i++)
	{	for (int j=0; j<v2.size(); j++)
		{	if (v1[i]==v2[j])
			{	in_both=true; break;
			}
		}
		if (!in_both)
			subtract_one_entry(v_out, v1[i]);
		in_both=false;
	}
	return v_out;
}

vector<int> MeshNamespace::merge_vectors(const vector<int> &v1, const vector<int> &v2)
{	vector<int> v_out;
	for (int i=0; i<v1.size(); i++)
		add_one_entry(v_out, v1[i]);
	for (int i=0; i<v2.size(); i++)
		add_one_entry(v_out, v2[i]);
	return v_out;
}
