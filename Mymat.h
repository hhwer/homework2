
#ifndef __MYMAT_H
#define __MYMAT_H

#include  <mpi.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

class Mymat
{

	public:
		Mymat();
		Mymat(int l,int m,int n);
		Mymat(int l,int m,int n,int num);
		~Mymat();
		
		void rank(int myid, int _size);
		
		void inposition();
//		void inposition(int myid, int size);		

		void outposition();
//		void outposition(int myid, int size);

		int* inbuff_x(int i);
		int* inbuff_y(int i);
		int* inbuff_z(int i);
		int* outbuff_x(int i);
		int* outbuff_y(int i);
		int* outbuff_z(int i);

		void myprint(void);




		std::vector<int> in_position;
		std::vector<int> xorder;
		std::vector<int> yorder;
		std::vector<int> zorder;
		std::vector<int> out_position_x;
		std::vector<int> out_position_y;
		std::vector<int> out_position_z;
		std::vector<int> ele;
	
	private:
		int size;
		int size_l;
		int size_m;
		int size_n;
		int myorder[3];
};

#endif

