#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>


class AdvectionElement
{
friend std::ostream& operator<<(std::ostream& os, const AdvectionElement& AE)
{	
	os<<AE.X[0]<<" "<<AE.X[1]<<" "<<AE.U[0]<<" "<<AE.U[1];
	return os;
}
public:
	// Pointer to the left neighbour
	AdvectionElement *Left_neighbour_pt;
	// Pointer to the right neighbour
	AdvectionElement *Right_neighbour_pt;
	// Storage for the coordinates
	std::vector<double> X;
	// Storage for the unknowns
	std::vector<double> U;
	// Constructor: initialise the vectors to hold two entries.
	AdvectionElement(double x_0, double x_1, double u_0, double u_1)
	{
	// Resize the vectors to hold two entries each
	// FILL THIS IN
	X.resize(2); U.resize(2);
	X[0]=x_0; X[1]=x_1;
	U[0]=u_0; U[1]=u_1;
	}
	//Returns the global coordinate x and unknown u values at the local coordinate s.
	double interpolated_x(double s) {return X[0]*(1/2*(1-s))+X[1]*(1/2*(1+s));}
	double interpolated_u(double s) {return U[0]*(1/2*(1-s))+U[1]*(1/2*(1+s));}
}; //End of the class definition

int main()
{
	/*

	//testing out pointers.		
	AdvectionElement zeroth(0,0,0,0), first(1,1,1,1),second(2,2,2,2), *another;			
	another=&zeroth;		//another is a pointer to zeroth
	zeroth.X[0]=600;
	std::cout<<zeroth<<" "<<*another<<std::endl; 	//outputs the value at the address pointed to by
						 	//another. Eventhough the zeroth.X[0] was
							// changed after the pointer was assigned,
							// the value at the address another is the
							// changed zeroth.X[0].

	std::cout<<zeroth.X[0]<<" "<<another->X[0]<<std::endl; 	//This is how you use member functions
								// of the object pointed to by a pointer
	*/

	/*

	//testing out pointers within a class definition
	AdvectionElement zeroth(0,0,0,0), first(1,1,1,1), second(2,2,2,2);
	first.Left_neighbour_pt=&zeroth; //first includes a pointer which points to zeroth
	AdvectionElement another=*first.Left_neighbour_pt; 	//another takes the value pointed to by 
								//the pointer
	zeroth.X[0]=15;
	std::cout<<zeroth.X[0]<<std::endl;			//obviously this is now 15
	std::cout<<first.Left_neighbour_pt->X[0]<<std::endl;	//this is 15 - the value at the address 
								//pointed to by f.Lnp
	std::cout<<(another.X[0])<<std::endl;			//this is still 0 - the same value as 									//zeroth had when another was assigned 									//the value at the pointer f.Lnp

	*/

	AdvectionElement zeroth(0,0,0,0), first(1,1,1,1), second(2,2,2,2);
	static const double pi=3.1415926535897932384626433832795;
	double domain_start=0.0, domain_end=2*pi;
	int no_elements=10;
	AdvectionElement element_array [3] {zeroth,first,second};

	for (int i=0; i<no_elements; i++)
	{
		
	}
	return 0;
}
