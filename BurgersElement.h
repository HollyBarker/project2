#ifndef BURGERSELEMENT_H
#define BURGERSELEMENT_H

#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<fstream>
#include "MMatrix.h"
#include "MVector.h"
#include "AdvectionElement.h"

class BurgersElement : public AdvectionElement
{
public:
	double flux(double u)
	{
		return (1.0/2.0)*u*u;
	}
	double h(double a, double b)
	{
		// Find max|f'(xi)|=max|xi| within the element by discretising the element into 10 points and choosing the largest
		double maxdflux=0 xi;
		for (int i=0;i<10;i++)
		{
			xi=a+i*((b-a)/10.0);	
			if (abs(xi)>maxdflux){maxdflux=xi;}
		}
		return (1.0/2.0)*(flux(a)+flux(b))-(1.0/2.0)*maxdflux*(b-a);
	}
};

#endif
