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
	// Pointer to the element's left neighbour
	AdvectionElement *Left_neighbour_pt;
	// Pointer to the element's right neighbour
	AdvectionElement *Right_neighbour_pt;
	// Vector to store x_0 and x_1
	std::vector<double> X;
	// Vector to store u_0 and u_1
	std::vector<double> U;
	// Default constructor
	AdvectionElement(){} 
	// Constructor to initialise element with x_0, x_1, u_0, u_1
	AdvectionElement(double x_0, double x_1, double u_0, double u_1)
	{
	// Provide a size for the vectors and fill the entries from the constructor
	X.resize(2); U.resize(2);
	X[0]=x_0; X[1]=x_1;
	U[0]=u_0; U[1]=u_1;
	}
	// Returns linear approximation for x(s) and u(s) for any point within the element
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
	double element_length=(domain_end-domain_start)/no_elements;
	//AdvectionElement element_array [3] ={zeroth,first,second};

	double x_0, x_1, u_0, u_1;
	std::vector<AdvectionElement> element_array(no_elements);
	for (int i=0; i<no_elements; i++)
	{
		x_0=domain_start+i*element_length;
		x_1=x_1+element_length;
		u_0=1.5*std::sin(x_0);
		u_1=1.5*std::sin(x_1);
		AdvectionElement AE(x_0,x_1,u_0,u_1);
		element_array[i]=AE;
	}
	element_array[0].Right_neighbour_pt=&element_array[1];
	element_array[no_elements-1].Left_neighbour_pt=&element_array[no_elements-2];
	for (int i=1; i<no_elements-1; i++)
	{
		element_array[i].Left_neighbour_pt=&element_array[i-1];
		element_array[i].Right_neighbour_pt=&element_array[i+1];
	}
	return 0;
}
