//

#include "Mymat.h"

int main(int argc, char** argv)
{
	int n,Max,myid,totalsize,size,ob;
	double err[2],total[2];
	int N=pow(2,2);
	ob = 0;
	if(argc>1)
	{
		N = atoi(argv[1]);
		if(argc>2)
		{
			ob = atoi(argv[2]);
		}
	}



	MPI::Init(argc, argv);
	myid = MPI::COMM_WORLD.Get_rank();
	totalsize = MPI::COMM_WORLD.Get_size();
	size = pow(totalsize+1.0, 1.0/3);
	n = N/size;
	Max = 1000;
	double mu=1.5;
	double n3 = pow(N,3);


	Mymat mat0(n,n,n,10);
	Mymat mat1(N,n,n/size);
	Mymat U(n,n,n);
	Mymat F(n,n,n);
	mat0.rank(myid,size);
	mat1.rank(myid,size);
	F.rank(myid,size);
	U.rank(myid,size);
	mat0.createtype(n);
	F.getF(N);

	mat0.outposition();
	mat1.inposition();
	mat0.createfactor(N,mu);

	
	for(int j=0;j<Max;j++)
    {	
		U = mat0;
		mat0.trans_z(mat1);

		mat0^=3;
		mat0 = mat0-F-(U*mu);

		//3d-fft
		mat0.trans_x(mat1);
		mat1.fft();
		mat0.retrans_x(mat1);
		mat0.trans_y(mat1);
		mat1.fft();
		mat0.retrans_y(mat1);
		mat0.trans_z(mat1);
		mat1.fft();

		mat0.retrans_z(mat1);
		mat0/=n3;
		mat0.dividefactor();
		mat0.trans_x(mat1);
		//ifft
		mat1.ifft();

		mat0.retrans_x(mat1);
		mat0.trans_y(mat1);
		//ifft
		mat1.ifft();

		mat0.retrans_y(mat1);
		mat0.trans_z(mat1);
		//ifft
		mat1.ifft();

		mat0.retrans_z(mat1);

		err[0] = (U-mat0).norm_inf();
		err[1] = U.norm_inf();
		MPI::COMM_WORLD.Allreduce(err, total, 2, MPI_DOUBLE, MPI_MAX);
		if(myid==0)
		{
			std::cout << "relative err= " << total[0]/(total[1]) << " ,err="					<< total[0]	<< std::endl;
		}
		
		if(total[0]/total[1]<1e-20 || total[0]==0)
		{
			break;
		}
  	
	}
  
//	if(myid == ob){
//
//		std::cout << "n=" << n << "N=" << N << std::endl;
//		std::ofstream fout{"rmat.data"};
//		for (auto i=0;i<n;i++){
//			for (auto j=0;j<n;j++){
//				for(auto k=0;k<n;k++){
//				fout << mat0.ele[i*n*n + j*n + k][0] << "+i" << mat1.ele[i*n*N + j *N +k][1] << "\t";
//				}
//				fout << std::endl;
//			}
//			fout << "\n";
//		}	
//    }
//	if(myid == ob){
//		std::cout << "n=" << n << "N=" << N << std::endl;
//		std::ofstream fout{"rmat.data"};
//		for (auto i=0;i<n;i++){
//			for (auto j=0;j<n;j++){
//				for(auto k=0;k<n;k++){
//				fout << mat0.factor[i*n*n + j*n + k] << "\t";
//				}
//				fout << std::endl;
//			}
//			fout << "\n";
//		}	
//    }
	mat0.typefree();
	MPI::Finalize();
}

