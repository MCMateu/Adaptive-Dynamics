#include<iostream>
#include<cmath>
#include<unistd.h>
#include<iomanip>
#include <stdio.h>
#include <math.h>  
#include <ctime>

using namespace std;

//He we are solving the coexistence condition
//See paper "The coevolution of predator-prey interactions: ESSS and Red Queen Dynamics"
// Paul Marrow, Richard Law and C. Cannings

double alpha(double s1)
{
	double c7,c8,c9;
	double aux;
	//double u;

	//u=0.001;

	c7=2.0;
	c8=8.0;
	c9=10.0;

	aux=(c7-c8*s1+c9*s1*s1);
	
	return aux;	


}

double beta(double s1, double s2)
{

	double delta1,delta2,c2,c3,c4,c5,c6;
	double aux;
	
	//double u;

	//u=0.001;

	c2=0.6;
	c3=0.5;
	c4=0.22;	
	c5=0.5;
	c6=0.25;

	delta1=(s1-c3)/c4;
	delta2=(s2-c5)/c6;

	aux=exp(-delta1*delta1+2.0*c2*delta1*delta2-delta2*delta2);

	return aux;

}

double gamma(double s1, double s2)
{

	double delta1,delta2;
	double c1,c2,c3,c4,c5,c6;
	double aux;
	//double u;

	//u=0.001;

	c1=0.2;
	c2=0.6;
	c3=0.5;
	c4=0.22;	
	c5=0.5;
	c6=0.25;
	
	delta1=(s1-c3)/c4;
	delta2=(s2-c5)/c6;

	aux=c1*exp(-delta1*delta1+2.0*c2*delta1*delta2-delta2*delta2);

	return aux;	

	

}


int main (void)
{

	//Files
	FILE* f3;
	f3=fopen("coexistencia.txt","w");

	FILE* f1;
	f1=fopen("phi1.txt","w");

	FILE* f2;
	f2=fopen("phi2.txt","w");

	FILE* f4;
	f4=fopen("phi1-part2.txt","w");

	
	
	//VARIABLE DECLARATION
	int i,j,k,l;
	double ft1,ft2,r1,r2;
	double h1,h2,s1,s2;
	double ovalo1,ovalo2;
	double x1,x2,hx;
	double phi1,phi2;
	double delta1,delta2;

	int contador,flag;
	double c1,c2,c3,c4,c5,c6,c7,c8,c9;
	//double u;

	//u=0.001;

	h1=0.00025;
	h2=0.00025;
	hx=0.0001;

	r1=0.5;
	r2=0.05;

	c1=0.2;
	c2=0.6;
	c3=0.5;
	c4=0.22;	
	c5=0.5;
	c6=0.25;
	c7=2.0;
	c8=8.0;
	c9=10.0;
	
	
		
	
	//COEXISTENCE REGION
	flag=0;
	contador=1;
	for(i=0;i<=4000;i++)
	{

		s1=h1*i;

		for(j=0;j<=4000;j++)
		{
			s2=h2*j;
			x1=r2/gamma(s1,s2);
			ovalo1=r1-alpha(s1)*x1;
			if(abs(ovalo1)<0.005)
			{
				ovalo2=-r2+gamma(s1,s2)*x1;
				if(abs(ovalo2)<0.0000005)
				{
					fprintf(f3,"%f %f \n",s1,s2);
					cout<<s1<<"  "<<s2<<endl;
						
				}
			}
		}	
		
	}
			
				
						
					
					
					
	
	

	//NULLCLINE LINES (NULL SELECTION)


		//SELECTION OF THE PRET
			
			//OUTSIDE THE OVAL
			for(j=0;j<=4000;j++)
			{
				s2=h2*j;
				fprintf(f4,"%f %f \n",c8/(2.0*c9),s2);

			}
			
			//INSIDE THE OVAL
			for(i=0;i<=4000;i++)
			{
				s1=h1*i;
				delta1=(s1-c3)/c4;
				
					for(j=495;j<=2920;j++)
					{

						s2=h2*j;
						delta2=(s2-c5)/c6;
	
						//POPULATIONS IN EQUILIBRIUM
						x1=r2/gamma(s1,s2);
						x2=(r1-alpha(s1)*x1)/beta(s1,s2);
										
						phi1=-x1*(-c8+2.0*c9*s1)-x2*beta(s1,s2)*(-2.0*delta1/c4+2.0*(c2/c4)*delta2);
						//cout<<phi1<<endl;
						if(abs(phi1)<0.005 && x1!=0 && x2!=0)
						{
							cout<<s1<<"\t"<<s2<<"\t"<<x1<<"\t"<<x2<<"\t"<<phi1<<endl;						
							fprintf(f1,"%f %f \n",s1,s2);
						}
											

								
					}

							
			}
		
	
		//SELECTION OF THE PREDATOR
		for(i=0;i<=400;i++)
		{
			s1=h1*i;
			delta1=(s1-c3)/c4;
			for(j=0;j<=400;j++)
			{

				s2=h2*j;				
	
				
				delta2=(s2-c5)/c6;

				phi2=(2.0*c2*delta1/c6-2.0*delta2/c6);

				if(abs(phi2)<0.0005)
				{		
					fprintf(f2,"%f %f \n",s1,s2);
				}
				
				
			}
		}

}
