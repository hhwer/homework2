//
#include "Mymat.h"

Mymat::Mymat(void)
{
	size_l = 0;
	size_m = 0;
	size_n = 0;
	status = 0;
}

Mymat::Mymat(int l, int m, int n )
{
	size_l = l;
	size_m = m;
	size_n = n;
	status = 0;
	ele = (fftw_complex*) malloc( sizeof(fftw_complex)*l*m*n );
	for(int i=0;i<l*m*n;i++)
	{	
		ele[i][0] = 0.0;
		ele[i][1] = 0.0;
	}
}


Mymat::Mymat(int l, int m, int n,int num)
{
	size_l = l;
	size_m = m;
	size_n = n;
	status = 0;
	double num0 = (double)num; 
	ele = (fftw_complex*) malloc( sizeof(fftw_complex)*l*m*n );
	for(int i=0;i<l*m*n;i++)
	{
		ele[i][0] = num0;
		ele[i][1] = 0.0;
	}
}


Mymat::Mymat(Mymat& mat1)
{
	size_l = mat1.size_l;
	size_m = mat1.size_m;
	size_n = mat1.size_n;
	status = 0;
	ele = (fftw_complex*) malloc( sizeof(fftw_complex)*size_l*size_m*size_n );
	for(int i=0;i<size_l*size_m*size_n;i++)
	{
		ele[i][0] = 0.0;
		ele[i][1] = 0.0;
	}
}
Mymat::~Mymat(void)
	{
		free(ele);
		if(status)
		{
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
	}

void Mymat::rank(int _myid, int _size)
{
	myid = _myid;
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
	status = 1;
	 byte_type = MPI::DOUBLE.Create_vector(1, 2, 2);
	byte_type.Commit();

	 tensor1_type = byte_type.Create_vector(n*n/size, n, n*size);
	tensor1_type.Commit();

	 xtensor0_type = byte_type.Create_vector(1, n*n*n/size, 0);
	xtensor0_type.Commit();

	 ycolumn0_type = byte_type.Create_vector(n, 1, n);
	ycolumn0_type.Commit();

	 ymatrix0_type = ycolumn0_type.Create_hvector(n, 1, 2*sizeof(double));
	ymatrix0_type.Commit();

	 ytensor0_type = ymatrix0_type.Create_vector(n/size, 1, 1);
	ytensor0_type.Commit();
	
	
	 zcolumn0_type = byte_type.Create_vector(n, 1, n*n);
	zcolumn0_type.Commit();

	 zmatrix0_type = zcolumn0_type.Create_hvector(n, 1, 											2*sizeof(double));
	zmatrix0_type.Commit();

	 ztensor0_type = zmatrix0_type.Create_hvector(n/size, 1, 										2*n*sizeof(double));
	ztensor0_type.Commit();


}

void Mymat::createfactor(int n,double mu)
{
	factor.resize(size_l*size_m*size_n);
	if(2*myorder[2]<size)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i),2) - mu;
				}
			}
		}
	}
	else if(2*myorder[2]>size)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i-n),2) - mu;
				}
			}
		}

	}
	else
	{	
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l/2;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i),2) - mu;
				}
				for(int i=size_l/2;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i-n),2) - mu;
				}
			}
		}
	}

	if(2*myorder[1]<size)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											-pow((myorder[1]*size_m+j),2);
				}
			}
		}
	}
	else if(2*myorder[1]>size)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											-pow((myorder[1]*size_m+j-n),2);
				}
			}
		}

	}
	else
	{	
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l/2;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											-pow((myorder[1]*size_m+j),2);
				}
				for(int i=size_l/2;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											-pow((myorder[1]*size_m+j-n),2);
				}
			}
		}
	}

	if(2*myorder[0]<size)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											-pow((myorder[0]*size_n+k),2);
				}
			}
		}
	}
	else if(2*myorder[0]>size)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											-pow((myorder[0]*size_n+k-n),2);
				}
			}
		}

	}
	else
	{	
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l/2;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											-pow((myorder[0]*size_n+k),2);
				}
				for(int i=size_l/2;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											-pow((myorder[0]*size_n+k-n),2);
				}
			}
		}
	}
}

void Mymat::dividefactor(void)
{
	for(int i=0;i<size_l*size_m*size_n;i++)
	{
		ele[i][0] /= factor[i];
		ele[i][1] /= factor[i];
	}
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
    			outbuff_x(i), 1, xtensor0_type, xorder[i], 99,								mat1.inbuff_x(i), 1, tensor1_type, mat1.xorder[i],  99);
    	}
}

void Mymat::trans_y(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
		MPI::COMM_WORLD.Sendrecv(
			outbuff_y(i), 1, ytensor0_type, yorder[i], 99,								mat1.inbuff_y(i), 1, tensor1_type, mat1.yorder[i],  99);
	}
}

void Mymat::trans_z(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
		MPI::COMM_WORLD.Sendrecv(
			outbuff_z(i), 1, ztensor0_type, zorder[i], 99,								mat1.inbuff_z(i), 1, tensor1_type, mat1.zorder[i],  99);
	}
}

void Mymat::retrans_x(Mymat &mat1)
{
	for(int i=0;i<size;i++)
    {
    	MPI::COMM_WORLD.Sendrecv(																mat1.inbuff_x(i), 1, tensor1_type, mat1.xorder[i],  99,						outbuff_x(i), 1, xtensor0_type, xorder[i], 99);
    }
}

void Mymat::retrans_y(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
    	MPI::COMM_WORLD.Sendrecv(																mat1.inbuff_y(i), 1, tensor1_type, mat1.yorder[i],  99,						outbuff_y(i), 1, ytensor0_type, yorder[i], 99);
	}
}

void Mymat::retrans_z(Mymat &mat1)
{
	for(int i=0;i<size;i++)
	{
    	MPI::COMM_WORLD.Sendrecv(																mat1.inbuff_z(i), 1, tensor1_type, mat1.zorder[i],  99,						outbuff_z(i), 1, ztensor0_type, zorder[i], 99);
	}
}


Mymat& Mymat::operator+=(const Mymat& mat1)
{
	for(int i=0;i<size_l*size_m*size_n;i++)
	{
		ele[i][0] += mat1.ele[i][0];
		ele[i][1] += mat1.ele[i][0];
	}
	return *this;
}

Mymat& Mymat::operator/=(double alpha)
{
	for(int i=0;i<size_l*size_m*size_n;i++)
	{
		ele[i][0] /= alpha;
		ele[i][1] /= alpha;
	}
	return *this;
}

Mymat& Mymat::operator=(const Mymat& mat1)
{
	for(int i=0;i<size_l*size_m*size_n;i++)
	{	
		ele[i][0] = mat1.ele[i][0];
		ele[i][1] = mat1.ele[i][1];
	}
	return *this;
}

Mymat& Mymat::operator^=(int p)
{
	int val;
	for(int i=0;i<size_l*size_m*size_n;i++)
	{
		val = ele[i][0]*ele[i][0] + ele[i][1]*ele[i][1];
		ele[i][0] *= val;
		ele[i][1] *= val;
	}
	return *this;
}

Mymat Mymat::operator-(const Mymat& mat1) const
{
	Mymat mat2(size_l,size_m,size_n);
	for(int i=0;i<size_l*size_m*size_n;i++)
	{
		mat2.ele[i][0] = ele[i][0] - mat1.ele[i][0];
		mat2.ele[i][1] = ele[i][1] - mat1.ele[i][0];
	}
	return mat2;
}


Mymat Mymat::operator*(double alpha) const
{	
	Mymat mat2(size_l,size_m,size_n);
	for(int i=0;i<size_l*size_m*size_n;i++)
	{
		mat2.ele[i][0] = ele[i][0] * alpha;
		mat2.ele[i][1] = ele[i][1] * alpha;
	}
	return mat2;
}



void Mymat::getF(int N)
{
	for(int k=0;k<size_n;k++)
	{
		double c1 = cos(((myorder[0]*size_n+k)/N-1)*M_PI);
		double s1 = sin(((myorder[0]*size_n+k)/N-1)*M_PI);
		for(int j=0;j<size_m;j++)
		{
			double c2 = cos(((myorder[1]*size_m+j)/N-1)*M_PI);
			double s2 = sin(((myorder[1]*size_m+j)/N-1)*M_PI);
			for(int i=0;i<size_l;i++)
			{
				ele[i][0] = s1+s2+sin(((myorder[0]*size_l+i)/N-1)*M_PI);
				ele[i][1] = c1+c2+cos(((myorder[0]*size_l+i)/N-1)*M_PI);
			}
		}
	}

}
