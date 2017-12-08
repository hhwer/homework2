//

#include "Mymat.h"

int main(int argc, char** argv)
{
	int n,Max,myid,totalsize,size;
	int N=pow(2,3);

	
	MPI::Init(argc, argv);
	myid = MPI::COMM_WORLD.Get_rank();
	totalsize = MPI::COMM_WORLD.Get_size();
	size = pow(totalsize+1.0, 1.0/3);
	n = N/size;
	Max = 1;
	double mu=1.5;
	double n3 = pow(N,3);

	int aaa = argc; 
	while(aaa==1)
	{	
	};

	Mymat mat0(n,n,n,myid);
	Mymat mat1(N,n,n/size);
	Mymat U(n,n,n);
	Mymat F(n,n,n);
	mat0.rank(myid,size);
	mat1.rank(myid,size);
	F.rank(myid,size);
	U.rank(myid,size);
	mat0.createtype(n);
//	mat1.createtype(n);
//	F.createtype(n);
//	U.createtype(n);
	F.getF(N);

	mat0.outposition();
	mat1.inposition();
	mat0.createfactor(N,mu);

	auto byte_type = MPI::DOUBLE.Create_vector(1, 2, 2);
	byte_type.Commit();

	for(int j=0;j<Max;j++)
    {	
		U = mat0;
		
		mat0^=3;
//		mat0 = mat0-F;
		mat0 = mat0-F-(U*mu);
		mat0.trans_x(mat1);
		//fft
		mat0.retrans_x(mat1);
		mat0.trans_y(mat1);
		//fft
		mat0.retrans_y(mat1);
		mat0.trans_z(mat1);
		//fft
		mat0.retrans_z(mat1);
		mat0/=n3;
		mat0.dividefactor();
		mat0.trans_x(mat1);
		//ifft
		mat0.retrans_x(mat1);
		mat0.trans_y(mat1);
		//ifft
		mat0.retrans_y(mat1);
		mat0.trans_z(mat1);
		//ifft
		mat0.retrans_z(mat1);
	}
	
	if(myid == 4)
	{
	std::ofstream fout{"rmat.data"};
	for (auto i=0;i<(n/size);i++){
		for (auto j=0;j<n;j++){
			for(auto k=0;k<N;k++){
			fout << mat1.ele[i*n*N + j*N + k][0] << "+i" << mat1.ele[i*n*N + j *N +k][1] << "\t";
			}
			fout << std::endl;
		}
		fout << "\n";
	}	
	MPI::Finalize();
	}
}

