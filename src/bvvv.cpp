#include "bvvv.h"

void bvvv :: init(int Num_ele, vector<int>& N_max)
{
	num_ele = Num_ele;
	n_max.resize(Num_ele);
	alpha.resize(Num_ele);
	beta.resize(Num_ele);
	gamma.resize(Num_ele);
	delta.resize(Num_ele);
	//
	n_b.resize(Num_ele);
	r_b.resize(Num_ele);
	l_b.resize(Num_ele);
	theta_b.resize(Num_ele);
	phi_b.resize(Num_ele);
	//
	n_max = N_max;
	//
	for(size_t t1=0; t1<Num_ele; t1++)
	{
		n_b[t1].resize(Num_ele);
		r_b[t1].resize(Num_ele);
		l_b[t1].resize(Num_ele);
		theta_b[t1].resize(N_max[t1]);
		phi_b[t1].resize(N_max[t1]);
	}
}
