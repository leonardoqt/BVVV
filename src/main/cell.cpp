#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include "cell.h"

void cell :: init(int Num_ele, vector<int> Num_atom, double Nei_r_max)
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
	num_nei.resize(num_atom,0);
	list_nei.resize(num_atom);
	ele_nei.resize(num_atom);
	pos_nei.resize(num_atom);
	nei_r_max = Nei_r_max;
}

void cell :: init(int Num_ele, double Nei_r_max)
{
	num_ele = Num_ele;
	num_atom_e.resize(Num_ele);
	nei_r_max = Nei_r_max;
}

void cell :: read_param(vec Param[3])
{
	param[0] = Param[0];
	param[1] = Param[1];
	param[2] = Param[2];
	vol = (param[0]^param[1])*param[2];
	inv_param[0] = (param[1]^param[2])/vol;
	inv_param[1] = (param[2]^param[0])/vol;
	inv_param[2] = (param[0]^param[1])/vol;
}

void cell :: read_atom_pos(int ind, int type, vec& pos)
{
	type_atom[ind] = type;
	pos_atom[ind] = pos;
}

void cell :: read_datafile(ifstream& in)
{
	string label_cell = "PRIMVEC";
	string label_pos  = "PRIMCOORD";
	string tmp;
	//
	while(getline(in,tmp))
		if(tmp.find(label_cell) != string::npos)
			break;
	in>>param[0]>>param[1]>>param[2];
	vol = (param[0]^param[1])*param[2];
	inv_param[0] = (param[1]^param[2])/vol;
	inv_param[1] = (param[2]^param[0])/vol;
	inv_param[2] = (param[0]^param[1])/vol;
	//
	while(getline(in,tmp))
		if(tmp.find(label_pos) != string::npos)
			break;
	in>>num_atom>>ene_dft;
	getline(in,tmp);
	//
	for(size_t t1=0; t1<num_ele; t1++)
		num_atom_e[t1] = 0;
	pos_atom.resize(num_atom);
	type_atom.resize(num_atom);
	num_nei.resize(num_atom,0);
	list_nei.resize(num_atom);
	ele_nei.resize(num_atom);
	pos_nei.resize(num_atom);
	for(size_t t1=0; t1<num_atom; t1++)
	{
		in>>type_atom[t1]>>pos_atom[t1];
		getline(in,tmp);
		num_atom_e[type_atom[t1]]++;
	}
}

// if cell can be very small, use this version
///*
void cell :: gen_nei_list()
{
	vec dr, dr_coef;
	double dot_p;
	for(size_t t1=0; t1<num_atom-1; t1++)
	for(size_t t2=t1+1; t2<num_atom; t2++)
	for(int nx=-2; nx<3; nx++)
	for(int ny=-2; ny<3; ny++)
	for(int nz=-2; nz<3; nz++)
	{
		dr = pos_atom[t2] - pos_atom[t1] + param[0]*nx+param[1]*ny+param[2]*nz;
		if(dr.norm() < nei_r_max)
		{
			num_nei[t1]++;
			list_nei[t1].push_back(t2);
			ele_nei[t1].push_back(type_atom[t2]);
			pos_nei[t1].push_back(dr);
			num_nei[t2]++;
			list_nei[t2].push_back(t1);
			ele_nei[t2].push_back(type_atom[t1]);
			pos_nei[t2].push_back(dr*(-1));
		}
	}
}
//*/

// if cell at least 2x2x2, use this simple version
/*
void cell :: gen_nei_list()
{
	vec dr, dr_coef;
	double dot_p;
	for(size_t t1=0; t1<num_atom-1; t1++)
	for(size_t t2=t1+1; t2<num_atom; t2++)
	{
		dr = pos_atom[t2] - pos_atom[t1];
		for(size_t t3=0; t3<3; t3++)
		{
			dot_p = dr*inv_param[t3];
			dot_p = dot_p-floor(dot_p+0.5);
			dr_coef[t3] = dot_p;
		}
		dr = param[0]*dr_coef[0]+param[1]*dr_coef[1]+param[2]*dr_coef[2];
		if(dr.norm() < nei_r_max)
		{
			num_nei[t1]++;
			list_nei[t1].push_back(t2);
			ele_nei[t1].push_back(type_atom[t2]);
			pos_nei[t1].push_back(dr);
			num_nei[t2]++;
			list_nei[t2].push_back(t1);
			ele_nei[t2].push_back(type_atom[t1]);
			pos_nei[t2].push_back(dr*(-1));
		}
	}
}
*/

void cell :: retrive_nei(int ind, int& type, vector<int>& type2, vector<vec>& pos)
{
	type = type_atom[ind];
	type2 = ele_nei[ind];
	pos = pos_nei[ind];
}

void cell :: assign_ene_dft(double Ene_dft)
{
	ene_dft = Ene_dft;
}

void cell :: get_ene_bvvv(bvvv& model)
{
	ene_bvvv = 0;
	for(size_t t1=0; t1<num_atom; t1++)
	{
		ene_bvvv += model.ene(type_atom[t1],ele_nei[t1],pos_nei[t1]);
		//cout<<t1<<'\t'<<ene_bvvv<<endl;
	}
}

void cell :: print()
{
	cout<<"Cell parameter     : "<<endl<<param[0]<<endl<<param[1]<<endl<<param[2]<<endl;
	cout<<"Volume             : "<<vol<<endl;
	cout<<"Energy from DFT    : "<<ene_dft<<endl;
	cout<<"Number of elements : "<<num_ele<<endl;
	cout<<"Number of atoms    : "<<num_atom<<endl;
	cout<<"Atom lists         : "<<endl;
	for(size_t t1=0; t1<num_atom; t1++)
		cout<<setw(4)<<t1+1<<" Type: "<<type_atom[t1]<<" Position: "<<pos_atom[t1]<<endl;
	cout<<"Neighbors          :"<<endl;
	for(size_t t1=0; t1<num_atom; t1++)
	{
		cout<<"    Atom "<<setw(4)<<t1+1<<", number of neighbors: "<<num_nei[t1]<<endl;
		for(size_t t2=0; t2<num_nei[t1]; t2++)
			cout<<"        "<<setw(4)<<list_nei[t1][t2]+1<<" ele "<<ele_nei[t1][t2]<<pos_nei[t1][t2]<<endl;
	}
}
