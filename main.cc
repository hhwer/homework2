//

#include "Mymat.h"

int main(int argc, char** argv)
{
	int n,myid,totalsize,size;
	int N=pow(2,3);


	MPI::Init(argc, argv);
	myid = MPI::COMM_WORLD.Get_rank();
	totalsize = MPI::COMM_WORLD.Get_size();
	size = pow(totalsize+1.0, 1.0/3);
	n = N/size;


	int aaa = 1;
	while(aaa==0)
	{	
	};

	Mymat mat0(n,n,n,myid);
	Mymat mat1(N,n,n/size);
	mat0.rank(myid,size);
	mat1.rank(myid,size);
	mat0.outposition();
	mat1.inposition();

//	auto byte_type = MPI::DOUBLE.Create_vector(1, 2, 2);
	auto byte_type = MPI::INT.Create_vector(1, 1, 1);
	byte_type.Commit();

	auto tensor1_type = byte_type.Create_vector(n*n/size, n, N);
	tensor1_type.Commit();

	auto xtensor0_type = byte_type.Create_vector(1, n*n*n/size, 0);
	xtensor0_type.Commit();

	auto ycolumn0_type = byte_type.Create_vector(n, 1, n);
	ycolumn0_type.Commit();

//	auto ymatrix0_type = ycolumn0_type.Create_hvector(n, 1, 2*sizeof(double));
	auto ymatrix0_type = ycolumn0_type.Create_hvector(n, 1, sizeof(int));
	ymatrix0_type.Commit();

	auto ytensor0_type = ymatrix0_type.Create_vector(n/size, 1, 1);
	ytensor0_type.Commit();
	
	
	auto zcolumn0_type = byte_type.Create_vector(n, 1, n*n);
	zcolumn0_type.Commit();

//	auto zmatrix0_type = zcolumn0_type.Create_hvector(n, 1, 											2*sizeof(double));
	auto zmatrix0_type = zcolumn0_type.Create_hvector(n, 1, 											sizeof(int));
	zmatrix0_type.Commit();

//	auto ztensor0_type = zmatrix0_type.Create_hvector(n/size, 1, 										2*n*sizeof(double));
	auto ztensor0_type = zmatrix0_type.Create_hvector(n/size, 1, 										n*sizeof(int));
	ztensor0_type.Commit();


	for(int i=0;i<size;i++)
	{
		MPI::COMM_WORLD.Sendrecv(
			mat0.outbuff_z(i), 1, ztensor0_type, mat0.zorder[i], 99,					mat1.inbuff_z(i), 1, tensor1_type, mat1.zorder[i],  99
		);
	}
	byte_type.Free();
	tensor1_type.Free();
	xtensor0_type.Free();
	ycolumn0_type.Free();
	ymatrix0_type.Free();
	ytensor0_type.Free();
	zcolumn0_type.Free();
	zmatrix0_type.Free();
	ztensor0_type.Free();
	
	if(myid == 5)
	{
	std::ofstream fout{"rmat.data"};
	for (auto i=0;i<(n/size);i++){
		for (auto j=0;j<n;j++){
			for(auto k=0;k<N;k++){
			fout << mat1.ele[i*n*N + j*N + k] << "\t";
			}
			fout << std::endl;
		}
		fout << "\n";
	}	
	MPI::Finalize();
	}


}

