#ifndef __BVVV__
#define __BVVV__

/*

Energy for BVVV

E_i = alpha_i*(n_i-n_max_i)^2 + beta_i*sin(pi*n_i)^2 + gamma_i*(|W_i| - delta_i*|sum_m f_m(n_i)*e_im|)^2

n_i      =  sum_j n_ij
n_ij     =  n_b_ij*e^(-(|r_ij|/|r_b_ij|)^l_b_ij)
W_i      =  sum_j n_ij*r_ij
f_m(n_i) =  tanh(4(n_i+0.5-m))
e_im     -> (theta_b_im, phi_b_im)

i,j are index of atoms, but different element(pairs) have the same coefficient

*/

#include <vector>

using namespace std;

class bvvv
{
private:
	int num_ele;
	vector <int> n_max;
	vector <double> alpha;
	vector <double> beta;
	vector <double> gamma;
	vector <double> delta;
	vector < vector <double> > n_b;
	vector < vector <double> > r_b;
	vector < vector <double> > l_b;
	vector < vector <double> > theta_b;
	vector < vector <double> > phi_b;
public:
	void init(int Num_ele, vector<int>& N_max);
};
#endif
