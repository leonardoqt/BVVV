#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <vec.h>
#include <mc.h>
#include "cell.h"
#include "bvvv.h"

using namespace std;

/*
	input structure:
	name of structure file
	number of structures
	rmax 
	number of elements
	max number of neighbors
	mc checkpoint
	mc period
*/

int main()
{
	// inputs
	ifstream in, xsf;
	// xsf related variables
	string name_xsf;
	int num_xsf;
	// cell related variables
	vector<cell> sys;
	int num_ele;
	double r_max;
	vector<int> max_nei;
	// bvvv related variables
	bvvv model;
	double err_tot;
	// mc related variables
	mc mc_control;
	thermo_profile temperature;
	int checkpoint;
	int period;
	int num_param_kind;
	vector <int> num_param_each;
	vector <double> lambda;
	vector<vector <double> > param;
	double T_max, T_min;
	// tmporary variables
	double ene_tmp;

	//====== inputs ======
	// read from input file
	in.open("bvvv.in");
	in>>name_xsf>>num_xsf>>r_max>>num_ele;
	max_nei.resize(num_ele);
	for(size_t t1=0; t1<num_ele; t1++)
		in>>max_nei[t1];
	in>>checkpoint>>period;
	in.close();

	// read from xsf file
	xsf.open(name_xsf);
	sys.resize(num_xsf);
	for(auto& m1 : sys)
	{
		m1.init(num_ele,r_max);
		m1.read_datafile(xsf);
		m1.gen_nei_list();
	}
	xsf.close();
	//====== internal initialization ======
	// initialize bvvv
	model.init(num_ele,max_nei);
	model.rand_self_assign();
	model.init_param_to_mc(param);
	model.send_param_to_mc(param);

	// initialize mc
	// get number of parameters and initial lambda
	num_param_kind = param.size();
	num_param_each.resize(num_param_kind);
	lambda.resize(num_param_kind);
	for(size_t t1=0; t1<num_param_kind; t1++)
	{
		num_param_each[t1] = param[t1].size();
		lambda[t1] = 0.2;
	}
	lambda[0] = 50;
	// get initial energy
	err_tot = 0;
	for(size_t t1=0; t1<num_xsf; t1++)
	{
		sys[t1].get_ene_bvvv(model);
		err_tot += (sys[t1].ene_bvvv-sys[t1].ene_dft)*(sys[t1].ene_bvvv-sys[t1].ene_dft);
	}
	err_tot/=num_xsf;
	// initialize
	mc_control.init(num_param_kind,checkpoint,num_param_each,lambda,param,err_tot);
	temperature.init(period,err_tot/10,0.001);

	// run mc
	for(size_t t1=0; t1<10*period; t1++)
	{
		mc_control.gen_param_kind(param);
		model.read_param_from_mc(param);
		err_tot = 0;
		for(auto& m2 : sys)
		{
			m2.get_ene_bvvv(model);
			err_tot += (m2.ene_bvvv-m2.ene_dft)*(m2.ene_bvvv-m2.ene_dft);
		}
		err_tot/=num_xsf;
		temperature.gen_T(t1,"mixed");
		mc_control.evaluate(temperature,err_tot);
		cout<<t1<<'\t'<<err_tot<<endl;
	}

	return 0;
}
