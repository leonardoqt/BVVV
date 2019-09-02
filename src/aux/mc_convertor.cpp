#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <vec.h>

using namespace std;

int main()
{
	ifstream log,xsf;
	ofstream out;
	string label_cell = "PRIMVEC";
	int num_ele = 2;
	int num_atom;
	vector<string> ele(num_ele);
	int num_xsf;
	string tmp;
	vec v_tmp;
	string ele_tmp;
	double ene_tmp;

	// could change to merge all files in sub dirs by using c++17 features
	xsf.open("save_trial.xsf");
	log.open("log.dat");
	// find element symbol
	log>>tmp;
	for(size_t t1=0; t1<num_ele; t1++)
		log>>ele[t1];
	log.clear(); log.seekg(ios::beg);
	// find number of structures
	num_xsf = 0;
	while(getline(xsf,tmp))
		if(tmp.find(label_cell) != string::npos)
			num_xsf++;
	xsf.clear(); xsf.seekg(ios::beg);

	// loop to store all structures
	out.open("example.xsf");
	getline(log,tmp);
	for(size_t t1=0; t1<num_xsf; t1++)
	{
		getline(xsf,tmp);
		out<<tmp<<endl;
		xsf>>v_tmp; out<<v_tmp<<endl;
		xsf>>v_tmp; out<<v_tmp<<endl;
		xsf>>v_tmp; out<<v_tmp<<endl;
		getline(xsf,tmp);
		getline(xsf,tmp);
		out<<tmp<<endl;
		xsf>>num_atom; getline(xsf,tmp);
		for(size_t t2=0; t2<num_ele+1; t2++)
			log>>tmp;
		log>>ene_tmp; getline(log,tmp);
		out<<num_atom<<' '<<ene_tmp<<endl;
		for(size_t t2=0; t2<num_atom; t2++)
		{
			xsf>>ele_tmp>>v_tmp; getline(xsf,tmp);
			out<<distance(ele.begin(),find(ele.begin(),ele.end(),ele_tmp))<<v_tmp<<endl;
		}
	}
}
