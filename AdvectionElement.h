#ifndef ADVECTIONELEMENT_H
#define ADVECTIONELEMENT_H

#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<fstream>
#include "MMatrix.h"
#include "MVector.h"


class AdvectionElement
{
friend std::ostream& operator<<(std::ostream& os, const AdvectionElement& AE)
{	
	os<<AE.X[0]<<" "<<AE.X[1]<<" "<<AE.U[0]<<" "<<AE.U[1];
	return os;
}
public:
	// Pointer to the element's left neighbour
	AdvectionElement *Left_neighbour_pt;
	// Pointer to the element's right neighbour
	AdvectionElement *Right_neighbour_pt;
	// Vector to store x_0 and x_1
	MVector X;
	// Vector to store u_0 and u_1
	MVector U;
	// Vector to store the timestep jump
	MVector timestepped_U;
	// Default constructor
	AdvectionElement(){} 
	// Constructor to initialise element with x_0, x_1, u_0, u_1
	AdvectionElement(double x_0, double x_1, double u_0, double u_1)
	{
	// Provide a size for the vectors and fill the entries from the constructor
	X.resize(2); U.resize(2); timestepped_U.resize(2);
	X[0]=x_0; X[1]=x_1;
	U[0]=u_0; U[1]=u_1;
	}
	// Returns linear approximation for x(s) and u(s) for any point within the element
	double interpolated_x(double s) {return X[0]*(1.0/2.0*(1.0-s))+X[1]*(1.0/2.0*(1.0+s));}
	double interpolated_u(double s) {return U[0]*(1.0/2.0*(1.0-s))+U[1]*(1.0/2.0*(1.0+s));}
	// Calculate the flux
	virtual double flux(double u)
	{
		double C=1;
		return C*u;
	}
	// Calculate the integral of the flux function over the element
	// using the two-point Gauss rule
	double integrate_flux()
	{
		double gauss_u1=interpolated_u(-sqrt(1.0/3.0));
		double gauss_u2=interpolated_u(sqrt(1.0/3.0)); 
		return (flux(gauss_u1)+flux(gauss_u2));
	}
	virtual double h(double a, double b)
	{
		//f=Cu so max(f')=C=1.
		double C=1;
		return (1.0/2.0)*(flux(a)+flux(b))-(1.0/2.0)*C*(b-a); 
	}
	void timestep(double dt)
	{
		//Construct MMatrix object for (M^e)^{-1}
		MMatrix invM(2,2,-2); 
		invM(0,0)=4; invM(1,1)=4; 
		invM=(1/(X[1]-X[0]))*invM; 
		
		//Construct MVector object for F^e
		MVector F(2);
		F[0]=-1; F[1]=1;
		F=(1.0/2.0)*integrate_flux()*F;

		//Construct MVector object for H^e
		MVector H(2);
		double U1em1=Left_neighbour_pt->U[1]; //u_1^{e-1}
		double U0ep1=Right_neighbour_pt->U[0];//u_0^{e+1}

		H[0]=h(U1em1,U[0]);
		H[1]=-h(U[1],U0ep1);

		timestepped_U=U+dt*(invM*(F+H));
	}
};

#endif
