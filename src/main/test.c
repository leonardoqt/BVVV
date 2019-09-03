#include <iostream>
#include <fstream>
#include <vector>
#include <vec.h>
#include "cell.h"
#include "bvvv.h"

using namespace std;

int main()
{
	/*
	ifstream in;
	cell sys[1011];
	in.open("example.xsf");
	for(size_t t1=0; t1<1011; t1++)
	{
		sys[t1].init(2,3.5);
		sys[t1].read_datafile(in);
		sys[t1].gen_nei_list();
	}
	sys[1010].print();
	return 0;
	*/
	//----------------------cell assignement test-------------------------
	///*
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
	//*/
	//---------------------BVVV assignment test----------------------
	///*
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
			c1[t2] = t1+1+1.0*t2;
			c2[t2] = t1+1+2.0*t2;
			c3[t2] = t1+1+3.0*t2;
		}
		d1.resize(n_max[t1]);
		d2.resize(n_max[t1]);
		for(size_t t2=0; t2<n_max[t1]; t2++)
		{
			d1[t2] = t1+0.4*t2;
			d2[t2] = t1+0.5*t2;
		}
		if(t1 == 1)
		{
			c1[2] = 1;
			c2[2] = 2/1.18232155679;
			c3[2] = 1;
			d1[0] = 0;
			d1[1] = 1;
			d1[2] = d1[3] = d1[4] = d1[5] = 0.5;
			d2[0] = d2[1] = d2[2] = 0;
			d2[3] = 0.25;
			d2[4] = 0.5;
			d2[5] = 0.75;
			model1.assign(1,1,0,0,1,1,c1,c2,c3,d1,d2);
		}
		else
			model1.assign(t1,t1,t1+1,t1+2,t1+3,t1+4,c1,c2,c3,d1,d2);
	}
	model1.print();
	//*/
	//-----------------------BVVV ene test-----------------------
	///*
	int e_t;
	vector<int> e_t2;
	vector<vec> p_t;
	vector< vector<double> > mc_tmp;
	sys1.retrive_nei(1,e_t,e_t2,p_t);
	cout<<model1.ene(e_t,e_t2,p_t)<<endl;
	model1.init_param_to_mc(mc_tmp);
	model1.send_param_to_mc(mc_tmp);
	model1.read_param_from_mc(mc_tmp);
	cout<<"---"<<endl;
	model1.print();
	cout<<model1.ene(e_t,e_t2,p_t)<<endl;
	sys1.get_ene_bvvv(model1);
	//*/
	return 0;
}
