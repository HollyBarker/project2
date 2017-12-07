#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<fstream>
#include "MMatrix.h"
#include "MVector.h"
#include "AdvectionElement.h"

static const double pi=3.141592653589793;

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

	std::ofstream dataFile;
	dataFile.open("initialElementData2.txt");
	if (!dataFile) return 1;
	

	//			Part 1 :: Problem set up

	// Set up the variables for the domain
	double domain_start=0.0, domain_end=2*pi;
	int no_elements=100;
	double element_length=(domain_end-domain_start)/no_elements;

	// Initialise the AdvectionElement inputs
	double x_0, x_1, u_0, u_1;
	// Initialise a vector of AdvectionElements
	std::vector<AdvectionElement> element_vector(no_elements);
	// Loop over the elements
	for (int i=0; i<no_elements; i++)
	{
		// Get the start and end coordinates for the element and implement the initial condition
		x_0=domain_start+i*element_length;
		x_1=x_0+element_length;
		if (x_0<1) u_0=1;
		else u_0=0;
		if (x_1<1) u_1=1;
		else u_1=0;

		// Use the coordinates and initial unknowns to make an AdvectionElement object
		AdvectionElement AE(x_0,x_1,u_0,u_1);
		// Fill in the vector of elements
		element_vector[i]=AE;
		// Write the initial condition for the central point of each element to a file
		dataFile.width(15); dataFile<<AE.interpolated_x(0);
		dataFile.width(15); dataFile<<AE.interpolated_u(0)<<std::endl;
	}

	// Assign pointers from the first element to the second element and final element
	element_vector[0].Right_neighbour_pt=&element_vector[1];
	element_vector[0].Left_neighbour_pt=&element_vector[no_elements-1];
	// Assign pointers from the final element to the penultimate element and first element
	element_vector[no_elements-1].Left_neighbour_pt=&element_vector[no_elements-2];
	element_vector[no_elements-1].Right_neighbour_pt=&element_vector[0];

	// Assign pointers from the remaining elements to their neighbours 
	for (int i=1; i<no_elements-1; i++)
	{
		element_vector[i].Left_neighbour_pt=&element_vector[i-1];
		element_vector[i].Right_neighbour_pt=&element_vector[i+1];
	}
	

	//			Part 2 :: Perform the timestep

	// Open a new .txt file to save the data
	std::ofstream dataFile2;
	dataFile2.open("timesteppedElementData2.txt");
	if (!dataFile2) return 1;

	// Set up variables to track the time 
	int no_timesteps=10;
	double time=0.0, timestep_length=0.01;

	// Loop over the timesteps
	for (int j=1; j<=no_timesteps;j++)
	{
		// Get the time at this timestep
		time=j*no_timesteps;
		// Loop over the elements
		//std::cout<<element_vector[30].U[0]<<std::endl;
		for (int i=0;i<no_elements;i++)
		{
			// Perform the timestep function on each AdvectionElement object
			element_vector[i].timestep(timestep_length);
		}
		//std::cout<<element_vector[30].U[0]<<std::endl;
		// Loop over the elements
		for (int i=0; i<no_elements;i++)
		{
			// Update all of the unknowns to their timestepped values AFTER the timestep has been calculated for each element
			element_vector[i].U = element_vector[i].timestepped_U;
		}
		//std::cout<<"updated"<<element_vector[30].U[0]<<std::endl;
	}
	
	// Loop over the elements
	for (int i=1;i<no_elements;i++)
	{
		// Write the AdvectionElement data to a file
		AdvectionElement AE=element_vector[i];
		dataFile2.width(15); dataFile2<<AE.interpolated_x(0);
		dataFile2.width(15); dataFile2<<AE.interpolated_u(0)<<std::endl;
	}
	return 0;
}
