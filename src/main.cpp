#include <iostream>
#include <fstream>
#include <vector>
#include <vec.h>
#include "cell.h"

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
	sys1.print();
	return 0;
}
