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
	ele.resize(l*m*n,0);
}


Mymat::Mymat(int l, int m, int n,int num)
{
	size_l = l;
	size_m = m;
	size_n = n;
	ele.resize(l*m*n,num);
}


Mymat::~Mymat(void){}

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



int* Mymat::inbuff_x(int i)
{
	int j =i;
	j = (j+myorder[2])%size;
	return &ele[in_position[j]];
}

int* Mymat::inbuff_y(int i)
{
	int j =i;
	j = (j+myorder[1])%size;
	return &ele[in_position[j]];
}


int* Mymat::inbuff_z(int i)
{
	int j =i;
	j = (j+myorder[0])%size;
	return &ele[in_position[j]];
}


int* Mymat::outbuff_x(int i)
{
	return &ele[out_position_x[i]];
}


int* Mymat::outbuff_y(int i)
{
	return &ele[out_position_y[i]];
}


int* Mymat::outbuff_z(int i)
{
	return &ele[out_position_z[i]];
}


void Mymat::myprint(void)
{
	std::cout << "size l,m,n=" << size << size_l << size_m << size_n <<std::endl;
	std::cout <<"inposition: ";
	for(int i=0;i<size;i++)
		std::cout << in_position[i] << " ";
	std::cout << std::endl;
}



