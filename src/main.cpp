#include <iostream>
#include <fstream>
#include <vector>
#include <vec.h>
#include "cell.h"
#include "bvvv.h"

using namespace std;

int main()
{
	ifstream in;
	vec param[3];
	int num_ele;
	vector<int> num_atom;
	double r_max;
	cell sys1;
	int ttype;
	vec tpos;
	
	in.open("in.dat");
	in>>param[0]>>param[1]>>param[2];
	in>>num_ele;
	num_atom.resize(num_ele);
	for(size_t t1=0; t1<num_ele; t1++)
		in>>num_atom[t1];
	in>>r_max;
	//
	sys1.init(num_ele,num_atom,r_max);
	//
	sys1.read_param(param);
	//
	for(size_t t1=0; t1<sys1.num_atom; t1++)
	{
		in>>ttype>>tpos;
		sys1.read_atom_pos(t1,ttype,tpos);
	}
	//
	sys1.gen_nei_list();
//	sys1.print();
	//--------------------------BVVV test-------------------------
	bvvv model1;
	vector<int> n_max(num_ele);
	vector<double> c1(num_ele),c2(num_ele),c3(num_ele),d1,d2;
	for(size_t t1=0; t1<num_ele; t1++)
		in>>n_max[t1];
	model1.init(num_ele,n_max);
	for(size_t t1=0; t1<num_ele; t1++)
	{
		for(size_t t2=0; t2<num_ele; t2++)
		{
			c1[t2] = t1+t2;
			c2[t2] = t1+2*t2;
			c3[t2] = t1+3*t2;
		}
		d1.resize(n_max[t1]);
		d2.resize(n_max[t1]);
		for(size_t t2=0; t2<n_max[t1]; t2++)
		{
			d1[t2] = t1+4*t2;
			d2[t2] = t1+5*t2;
		}
		model1.assign(t1,t1+1,t1+2,t1+3,t1+4,c1,c2,c3,d1,d2);
	}
	model1.print();
	return 0;
}
