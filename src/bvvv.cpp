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

void bvvv :: assign(int Ele, double Alpha, double Beta, double Gamma, double Delta, vector<double>& N_b, vector<double>& R_b, vector<double>& L_b, vector<double>& Theta_b, vector<double>& Phi_b)
{
	if (Ele >= num_ele)
	{
		cout<<"Error in assignment, invalid element index "<<Ele<<endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		alpha[Ele] = Alpha;
		beta[Ele]  = Beta;
		gamma[Ele] = Gamma;
		delta[Ele] = Delta;
		if (N_b.size()!=num_ele || R_b.size()!=num_ele || L_b.size()!=num_ele || Theta_b.size()!=n_max[Ele] || Phi_b.size()!=n_max[Ele])
		{
			cout<<"Error in assignemnt, wrong length of paired parameters"<<endl;
			exit(EXIT_FAILURE);
		}
		else
		{
			n_b[Ele] = N_b;
			r_b[Ele] = R_b;
			l_b[Ele] = L_b;
			theta_b[Ele] = Theta_b;
			phi_b[Ele] = Phi_b;
		}
	}
}


void bvvv :: print()
{
	cout<<"Number of elements        : "<<num_ele<<endl;
	cout<<"Maximum N of each element : "<<endl;
	for(size_t t1=0; t1<num_ele; t1++) cout<<'\t'<<n_max[t1]; cout<<endl;
	cout<<"Alpha   :"<<endl;
	for(size_t t1=0; t1<num_ele; t1++) cout<<'\t'<<alpha[t1]; cout<<endl;
	cout<<"Beta    :"<<endl;
	for(size_t t1=0; t1<num_ele; t1++) cout<<'\t'<<beta[t1]; cout<<endl;
	cout<<"Gamma   :"<<endl;
	for(size_t t1=0; t1<num_ele; t1++) cout<<'\t'<<gamma[t1]; cout<<endl;
	cout<<"Delta   :"<<endl;
	for(size_t t1=0; t1<num_ele; t1++) cout<<'\t'<<delta[t1]; cout<<endl;
	cout<<"N in BV :"<<endl;
	for(size_t t1=0; t1<num_ele; t1++)
	{
		for(size_t t2=0; t2<num_ele; t2++)
			cout<<'\t'<<n_b[t1][t2];
		cout<<endl;
	}
	cout<<"R in BV :"<<endl;
	for(size_t t1=0; t1<num_ele; t1++)
	{
		for(size_t t2=0; t2<num_ele; t2++)
			cout<<'\t'<<r_b[t1][t2];
		cout<<endl;
	}
	cout<<"L in BV :"<<endl;
	for(size_t t1=0; t1<num_ele; t1++)
	{
		for(size_t t2=0; t2<num_ele; t2++)
			cout<<'\t'<<l_b[t1][t2];
		cout<<endl;
	}
	cout<<"Theta in BVVV :"<<endl;
	for(size_t t1=0; t1<num_ele; t1++)
	{
		for(size_t t2=0; t2<n_max[t1]; t2++)
			cout<<'\t'<<theta_b[t1][t2];
		cout<<endl;
	}
	cout<<"Phi in BVVV :"<<endl;
	for(size_t t1=0; t1<num_ele; t1++)
	{
		for(size_t t2=0; t2<n_max[t1]; t2++)
			cout<<'\t'<<phi_b[t1][t2];
		cout<<endl;
	}
}
