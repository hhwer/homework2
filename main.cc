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
	Max = 100;

	int aaa = 1;
	while(aaa==0)
	{	
	};

	Mymat mat0(n,n,n,myid);
	Mymat mat1(N,n,n/size);
	mat0.rank(myid,size);
	mat1.rank(myid,size);
	mat0.outposition();
	mat0.createtype(n);
	mat1.inposition();

//	auto byte_type = MPI::DOUBLE.Create_vector(1, 2, 2);
//	byte_type.Commit();
//
//	auto tensor1_type = byte_type.Create_vector(n*n/size, n, N);
//	tensor1_type.Commit();
//
//	auto xtensor0_type = byte_type.Create_vector(1, n*n*n/size, 0);
//	xtensor0_type.Commit();
//
//	auto ycolumn0_type = byte_type.Create_vector(n, 1, n);
//	ycolumn0_type.Commit();
//
//	auto ymatrix0_type = ycolumn0_type.Create_hvector(n, 1, 2*sizeof(double));
//	ymatrix0_type.Commit();
//
//	auto ytensor0_type = ymatrix0_type.Create_vector(n/size, 1, 1);
//	ytensor0_type.Commit();
//	
//	
//	auto zcolumn0_type = byte_type.Create_vector(n, 1, n*n);
//	zcolumn0_type.Commit();
//
//	auto zmatrix0_type = zcolumn0_type.Create_hvector(n, 1, 											2*sizeof(double));
//	zmatrix0_type.Commit();
//
//	auto ztensor0_type = zmatrix0_type.Create_hvector(n/size, 1, 										2*n*sizeof(double));
//	ztensor0_type.Commit();

	for(int j=0;j<Max;j++)
    {
//		for(int i=0;i<size;i++)
//    	{
//    		MPI::COMM_WORLD.Sendrecv(
//    			mat0.outbuff_x(i), 1, mat0.xtensor0_type, mat0.xorder[i], 99,					mat1.inbuff_x(i), 1, mat0.tensor1_type, mat1.xorder[i],  99
//    		);
//    	}
		mat0.trans_x(mat1);
		mat0.trans_y(mat1);
		mat0.trans_z(mat1);

//		for(int i=0;i<size;i++)
//    	{
//    		MPI::COMM_WORLD.Sendrecv(
//    			mat0.outbuff_y(i), 1, mat0.ytensor0_type, mat0.yorder[i], 99,					mat1.inbuff_y(i), 1, mat0.tensor1_type, mat1.yorder[i],  99
//    		);
//    	}
//
//		for(int i=0;i<size;i++)
//    	{
//    		MPI::COMM_WORLD.Sendrecv(
//    			mat0.outbuff_z(i), 1, mat0.ztensor0_type, mat0.zorder[i], 99,					mat1.inbuff_z(i), 1, mat0.tensor1_type, mat1.zorder[i],  99
//    		);
//    	}
	}
//	byte_type.Free();
//	tensor1_type.Free();
//	xtensor0_type.Free();
//	ycolumn0_type.Free();
//	ymatrix0_type.Free();
//	ytensor0_type.Free();
//	zcolumn0_type.Free();
//	zmatrix0_type.Free();
//	ztensor0_type.Free();
	
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

