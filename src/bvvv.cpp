#include <cmath>
#include "bvvv.h"

void bvvv :: init(int Num_ele, vector<int>& N_max)
{
	if(Num_ele != N_max.size())
	{
		cout<<"Error in bvvv initialization, wrong number of n_max for each element"<<endl;
		exit(EXIT_FAILURE);
	}
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
	// guard
	if (Ele >= num_ele)
	{
		cout<<"Error in assignment, invalid element index "<<Ele<<endl;
		exit(EXIT_FAILURE);
	}
	//
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

double bvvv :: ene(int ele0, vector<int>& ele1, vector<vec>& bond)
{
	if(ele1.size() != bond.size())
	{
		cout<<"Error in calculating energy, number of atoms not consistent"<<endl;
		exit(EXIT_FAILURE);
	}
	int num = ele1.size();
	vector <double> n_j(num);
	vec ww,ww0,tmp;
	double sum_n, theta, phi;
	double etot=0;
	// n_j, n, ww
	ww = 0.0;
	sum_n = 0;
	for(size_t t1=0; t1<num; t1++)
	{
		n_j[t1] = n_b[ele0][ele1[t1]]*exp(1-pow(bond[t1].norm()/fabs(r_b[ele0][ele1[t1]]),fabs(l_b[ele0][ele1[t1]])));
		sum_n += n_j[t1];
		ww = ww + bond[t1]*n_j[t1];
		//cout<<bond[t1]<<'\t'<<ele1[t1]<<'\t'<<'\t'<<'\t'<<n_j[t1]<<endl;
	}
	// ww0
	ww0 = 0.0;
	for(size_t t1=0; t1<n_max[ele0]; t1++)
	{
		theta = M_PI*theta_b[ele0][t1];
		phi = 2*M_PI*phi_b[ele0][t1];
		tmp[0] = sin(theta)*cos(phi);
		tmp[1] = sin(theta)*sin(phi);
		tmp[2] = cos(theta);
		ww0 = ww0 + tmp*((1+tanh(4*(sum_n-0.5-t1)))/2);
	}
	//cout<<"W"<<ww<<"\tW0"<<ww0<<endl;
	// ene
	etot = alpha[ele0]*pow(sum_n-n_max[ele0],2)+beta[ele0]*pow(sin(M_PI*sum_n),2)+gamma[ele0]*pow(ww.norm() - delta[ele0]*ww0.norm(),2);
	return etot;
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
