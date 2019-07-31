#ifndef __CELL__
#define __CELL__

#include <vec.h>
#include <vector>

class cell
{
private:
	vec param[3];
	int num_ele;
	int num_atom;
	vector<int> num_atom_e;
	vector<vec> pos_atom;
	vector<int> type_atom;
	vector<int> num_nei;
	vector< vector<int> > list_nei;
	vector< vector<vec> > corr_nei;
	double nei_r_max;
public:
	// interface related (will change for real code)
	void init(int Num_ele, vector<int> Num_atom, double Nei_r_max);
	void read_param(vec Param[3]);
	void read_atom_pos(int ind, int type, vec& pos);
};

#endif
