#include "cell.h"

void cell :: init(int Num_ele, vec<int> Num_atom, double Nei_r_max)
{
	num_ele = Num_ele;
	num_atom_e.resize(Num_ele);
	num_atom = 0;
	for(size_t t1=0; t1<Num_ele; t1++)
	{
		num_atom_e[t1] = Num_atom[t1];
		num_atom += Num_atom[t1];
	}
	pos_atom.resize(num_atom);
	type_atom.resize(num_atom);
	num_nei.resize(num_atom);
	list_nei.resize(num_atom);
	corr_nei.resize(num_atom);
	nei_r_max = Nei_r_max;
}

void cell :: read_param(vec Param[3])
{
	param[0] = Param[0];
	param[1] = Param[1];
	param[2] = Param[2];
}

void cell :: read_atom_pos(int ind, int type, vec& pos)
{
	
}
