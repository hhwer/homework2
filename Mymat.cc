//
#include "Mymat.h"

Mymat::Mymat(void)
{
	size_l = 0;
	size_m = 0;
	size_n = 0;
}

Mymat::Mymat(int l, int m, int n )
{
	size_l = l;
	size_m = m;
	size_n = n;
	ele = (fftw_complex*) malloc( sizeof(fftw_complex)*l*m*n );
	for(int i=0;i<l*m*n;i++)
	{	
		ele[i][0] = 0.0;
		ele[i][1] = 0.0;
	}
//	ele = fftw_alloc_complex(l*m*n);
}


Mymat::Mymat(int l, int m, int n,int num)
{
	size_l = l;
	size_m = m;
	size_n = n;
	double num0 = (double)num; 
	ele = (fftw_complex*) malloc( sizeof(fftw_complex)*l*m*n );
//	ele = fftw_alloc_complex(l*m*n);
	for(int i=0;i<l*m*n;i++)
	{
		ele[i][0] = num0;
		ele[i][1] = 0.0;
	}
}


Mymat::~Mymat(void)
	{
		free(ele);
		byte_type.Free();
		tensor1_type.Free();
		xtensor0_type.Free();
		ycolumn0_type.Free();
		ymatrix0_type.Free();
		ytensor0_type.Free();
		zcolumn0_type.Free();
		zmatrix0_type.Free();
		ztensor0_type.Free();
	}

void Mymat::rank(int myid, int _size)
{
	size = _size;
	xorder.resize(size);
	yorder.resize(size);
	zorder.resize(size);

	myorder[2] = myid%size;
	myorder[0] = myid/size/size;
	myorder[1] = myid/size%size;
	for(int i=0;i<size;i++)
	{
		xorder[i] = (myorder[2]+i)%size + myorder[1] * size +												 myorder[0] * size * size;
		yorder[i] = myorder[2] + (myorder[1]+i)%size * size +												 myorder[0] * size * size;
		zorder[i] = myorder[2] + myorder[1] * size +												 (myorder[0]+i)%size * size * size;
	}
	
}



void Mymat::inposition(void)
{
	in_position.resize(size);
	for(int i=0;i<size;i++)
	{
		in_position[i] = i*size_l/size; 
	}
}

void Mymat::outposition(void)
{
	out_position_x.resize(size);
	out_position_y.resize(size);
	out_position_z.resize(size);
	
	for(int i=0;i<size;i++)
	{
		out_position_x[i] = i*size_l*size_m*size_n/size; 
		out_position_y[i] = i*size_l*size_m*size_n/size; 
		out_position_z[i] = i*size_l; 
	}
}

void Mymat::createtype(int n)

{
	auto byte_type = MPI::DOUBLE.Create_vector(1, 2, 2);
	byte_type.Commit();

	auto tensor1_type = byte_type.Create_vector(n*n/size, n, n*size);
	tensor1_type.Commit();

	auto xtensor0_type = byte_type.Create_vector(1, n*n*n/size, 0);
	xtensor0_type.Commit();

	auto ycolumn0_type = byte_type.Create_vector(n, 1, n);
	ycolumn0_type.Commit();

	auto ymatrix0_type = ycolumn0_type.Create_hvector(n, 1, 2*sizeof(double));
	ymatrix0_type.Commit();

	auto ytensor0_type = ymatrix0_type.Create_vector(n/size, 1, 1);
	ytensor0_type.Commit();
	
	
	auto zcolumn0_type = byte_type.Create_vector(n, 1, n*n);
	zcolumn0_type.Commit();

	auto zmatrix0_type = zcolumn0_type.Create_hvector(n, 1, 											2*sizeof(double));
	zmatrix0_type.Commit();

	auto ztensor0_type = zmatrix0_type.Create_hvector(n/size, 1, 										2*n*sizeof(double));
	ztensor0_type.Commit();


}

double* Mymat::inbuff_x(int i)
{
	int j =i;
	j = (j+myorder[2])%size;
	return &(ele[in_position[j]][0]);
}

double* Mymat::inbuff_y(int i)
{
	int j =i;
	j = (j+myorder[1])%size;
	return &(ele[in_position[j]][0]);
}


double* Mymat::inbuff_z(int i)
{
	int j =i;
	j = (j+myorder[0])%size;
	return &(ele[in_position[j]][0]);
}


double* Mymat::outbuff_x(int i)
{
	return &(ele[out_position_x[i]][0]);
}


double* Mymat::outbuff_y(int i)
{
	return &(ele[out_position_y[i]][0]);
}


double* Mymat::outbuff_z(int i)
{
	return &(ele[out_position_z[i]][0]);
}


void Mymat::trans_x(Mymat &mat1)
{
	for(int i=0;i<size;i++)
    {
    	MPI::COMM_WORLD.Sendrecv(
    			outbuff_x(i), 1, xtensor0_type, xorder[i], 99,					mat1.inbuff_x(i), 1, tensor1_type, mat1.xorder[i],  99
    		);
    	}
}

void Mymat::trans_y(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
		MPI::COMM_WORLD.Sendrecv(
			outbuff_y(i), 1, ytensor0_type, yorder[i], 99,					mat1.inbuff_y(i), 1, tensor1_type, mat1.yorder[i],  99
		);
	}
}

void Mymat::trans_z(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
		MPI::COMM_WORLD.Sendrecv(
			outbuff_z(i), 1, ztensor0_type, zorder[i], 99,					mat1.inbuff_z(i), 1, tensor1_type, mat1.zorder[i],  99
		);
	}
}





