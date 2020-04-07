#include<iostream>
#include<cmath>
#include<unistd.h>
#include<iomanip>
#include <stdio.h>
#include <math.h>  
#include <ctime>

using namespace std;

//The equations solved here are the ones presented in the paper:
//The dynamical theory of coevolution: a derivation from stochastic ecological processes, Ulf Dieckmann & Richard Laq
//You can change the parameter set, some ideas are provided in the paper.

//EXTERNAL FUNCTIONS ARE DEFFINED ACCORDING WITH PAPER

//Intraspecific competence for the preys
double alpha(double s1)
{
	double c7,c8,c9;
	double aux;
	double u;

	u=0.001;

	c7=2.0;
	c8=8.0;
	c9=10.0;

	aux=(c7-c8*s1+c9*s1*s1);
	
	return aux*u;	


}
//Predation term for the preys
double beta(double s1, double s2)
{

	double delta1,delta2,c2,c3,c4,c5,c6;
	double aux;
	
	double u;

	u=0.001;

	c2=0.6;
	c3=0.5;
	c4=0.22;	
	c5=0.5;
	c6=0.25;

	delta1=(s1-c3)/c4;
	delta2=(s2-c5)/c6;

	aux=exp(-delta1*delta1+2.0*c2*delta1*delta2-delta2*delta2);

	return aux*u;

}
//Predation term for the predator (benefit for the predators)
double gamma(double s1, double s2)
{

	double delta1,delta2;
	double c1,c2,c3,c4,c5,c6;
	double aux;
	double u;

	u=0.001;

	c1=0.2;
	c2=0.6;
	c3=0.5;
	c4=0.22;	
	c5=0.5;
	c6=0.25;
	
	delta1=(s1-c3)/c4;
	delta2=(s2-c5)/c6;

	aux=c1*exp(-delta1*delta1+2.0*c2*delta1*delta2-delta2*delta2);

	return aux*u;	

	

}


double M1(double s1 , double As1 )
{
	double pi,sigma1;
	double aux;
	
	pi=3.14159265359;
	sigma1=5.0E-3;

	aux=(1.0/(sqrt(2.0*pi)*sigma1))*exp(-0.5*(As1*As1)/(sigma1*sigma1));

	return aux;

} 

double M2(double s2 , double As2 )
{
	double pi,sigma2;
	double aux;

	
	pi=3.14159265359;
	sigma2=5.0E-3;

	aux=(1.0/(sqrt(2.0*pi)*sigma2))*exp(-0.5*(As2*As2)/(sigma2*sigma2));

	return aux;

} 


int main (void)
{

	FILE* f1;
	f1=fopen("trayectoria4.txt","w");
	

	//VARIABLE DECLARATION
	int i,j,k,l;
	double Ns,Nt;
	double ft1,ft2,r1,r2;
	double h1,h2,s1,s2;
	double hx,t;
	double phi1,phi2;
	double delta1,delta2;
	double s1i,s1f,s2i,s2f;
	double hs,ht;	

	int contador,flag;
	double c1,c2,c3,c4,c5,c6,c7,c8,c9;
	double mu1,mu2;
	double n1,n2;
	double sigma1,sigma2;
	double courant;

	//PARAMETERS
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

	mu1=1.0E-4;
	mu2=1.0E-3;

	sigma1=5.0E-3;
	sigma2=5.0E-3;


	//INTEGRATION PARAMTERS
	Ns=10000;
	Nt=1E10;
	
	ht=1.0E6/(1.0*Nt);
	hs=(1.0-0.0)/(1.0*Ns);
	
	courant=ht/(hs*hs);

	cout<<"Número de courant:   "<<courant<<endl;

	//PVI
	s1i=0.35;
	s2i=0.16;
	
	fprintf(f1,"%f %f \n",s1i,s2i);

	contador=0;
	flag=1E6;
	for(i=1;i<=Nt;i++)
	{
		

		//EQUILIBRIUM OF THE NUMER OF INDIVIDUAL OF EACH SPECIES
		n1=r2/gamma(s1i,s2i);
		n2=(r1-alpha(s1i)*n1)/beta(s1i,s2i);
		
		//TRAITS DERIVATIVES (CANONICAL EQUATIONS)	

		
		s1f=s1i+ht*0.5*mu1*sigma1*sigma1*n1*(r1-alpha(s1i+hs)*n1-beta(s1i+hs,s2i)*n2)/(hs);

		s2f=s2i+ht*0.5*mu2*sigma2*sigma2*n2*(-r2+n1*gamma(s1i,s2i+hs))/(hs);


		//RESET VALUES
		s1i=s1f;
		s2i=s2f;	

				contador++;
				if(contador==flag)
				{
					cout<<setprecision(15)<<s1f<<"\t"<<s2f<<endl;
					
					fprintf(f1,"%f %f \n",s1i,s2i);

			
					contador=0;
				}
	}	



}
