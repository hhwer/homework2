
#include "Mymat.h"


/* --------------------------------------------------------------------------*/
/**
* @brief 从\Delta U-\mu U的fourier系数得到U的傅里叶系数所需的除数
*
* @param N  单个方向的规模
* @param mu	
*/
/* ----------------------------------------------------------------------------*/
void Mymat::createfactor(int N,double mu)
{
	factor.resize(size_l*size_m*size_n);
	if(2*myorder[2]<size-1)
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
	else if(2*myorder[2]>size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i-N),2) - mu;
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
					factor[i+j*size_l+k*size_l*size_m] = 											-pow((myorder[2]*size_l+i-N),2) - mu;
				}
			}
		}
	 }
	

	if(2*myorder[1]<size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[1]*size_m+j),2);
				}
			}
		}
	}
	else if(2*myorder[1]>size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[1]*size_m+j-N),2);
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
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[1]*size_m+j),2);
				}
				for(int i=size_l/2;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[1]*size_m+j-N),2);
				}
			}
		}
	}

	if(2*myorder[0]<size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[0]*size_n+k),2);
				}
			}
		}
	}
	else if(2*myorder[0]>size-1)
	{
		for(int k=0;k<size_n;k++)
		{
			for(int j=0;j<size_m;j++)
			{
				for(int i=0;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[0]*size_n+k-N),2);
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
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[0]*size_n+k),2);
				}
				for(int i=size_l/2;i<size_l;i++)
				{
					factor[i+j*size_l+k*size_l*size_m] -= 											pow((myorder[0]*size_n+k-N),2);
				}
			}
		}
	}
}

/* --------------------------------------------------------------------------*/
/**
* @brief 从\Delta U -\mu U的fourier得到U的fourier
*/
/* ----------------------------------------------------------------------------*/
void Mymat::dividefactor(void)
{
	for(int i=0;i<size_l*size_m*size_n;i++)
	{
		ele[i][0] /= factor[i];
		ele[i][1] /= factor[i];
	}
}


/* --------------------------------------------------------------------------*/
/**
* @brief 生成右端项
*
* @param N 问题单个方向上的规模
*/
/* ----------------------------------------------------------------------------*/
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
				ele[i+j*size_l+k*size_l*size_m][0] = 												s1*s2*sin(((myorder[2]*size_l+i)/N-1)*M_PI);
				ele[i+j*size_l+k*size_l*size_m][1] = 												c1*c2*cos(((myorder[2]*size_l+i)/N-1)*M_PI);
			}
		}
	}

}

void Mymat::getvalue(void)
{
	for(int k=0;k<size_n;k++)
	{
		for(int j=0;j<size_m;j++)
		{
			for(int i=0;i<size_l;i++)
			{
			ele[i+j*size_l+k*size_l*size_m][0] = i*myid;
			}
		}
	}

}


