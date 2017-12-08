
#include "Mymat.h"


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


