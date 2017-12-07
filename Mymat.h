
#ifndef __MYMAT_H
#define __MYMAT_H

#include <fftw3.h>
#include <mpi.h>
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
		void outposition();
		void createtype(int n);

		void trans_x(Mymat &mat1);
		void trans_y(Mymat &mat1);
		void trans_z(Mymat &mat1);

		fftw_complex* ele;

	protected:
		double* inbuff_x(int i);
		double* inbuff_y(int i);
		double* inbuff_z(int i);
		double* outbuff_x(int i);
		double* outbuff_y(int i);
		double* outbuff_z(int i);




		void myprint(void);


		

	
	private:
		int size;
		int size_l;
		int size_m;
		int size_n;
		int myorder[3];
		std::vector<int> in_position;
		std::vector<int> out_position_x;
		std::vector<int> out_position_y;
		std::vector<int> out_position_z;
		std::vector<int> xorder;
		std::vector<int> yorder;
		std::vector<int> zorder;

		MPI::Datatype byte_type;
		MPI::Datatype tensor1_type;
		MPI::Datatype xtensor0_type;
		MPI::Datatype ycolumn0_type;
		MPI::Datatype ymatrix0_type;
		MPI::Datatype ytensor0_type;
		MPI::Datatype zcolumn0_type;
		MPI::Datatype zmatrix0_type;
		MPI::Datatype ztensor0_type;
};

#endif

