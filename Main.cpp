//#define _CRT_SECURE_NO_DEPRECATE
#include "Lattice.h"
#include<iostream>
#include<fstream>
#include<time.h>
using namespace std;

int main(){
	int M = 3;
	int y;		// No use. Just to stop the execution at any point using cin>>y
	long seed;
	long seed2;
	int latticeLength = 35;			//Total no. of particles in a row
	int latticeBreadth = 35;			//Total no. of particles in a column     B=2*L-1 for constraint hexagon
	int NoOfEvolutions = 1;			//Total no. of repeated simulations
	int totalIterations = 100000;		//Total no. of iterations in each simulation
	int QuiverplotInterval = 100;		// according to totalIterations and h
	int forceHistInterval = 500;
	int StartRowNo = 1;					//for Histogram plot of force
	int RowNoInterval = 3;				//for Histogram plot of force
	int RowNo = StartRowNo;				//for Histogram plot of force
	int count_RowNo = 1;				//for Histogram plot of force
	int falseInput;						// No use of this in system/evolution; only used in continuing the execution of code
	//	Particle Particle_checkEvolution;
//	int particleNo;	
//	int count = 0;				// used in checkEvolution f:n

	Lattice hex = Lattice(latticeLength, latticeBreadth, 1.0);  // last no. is lattice length. MAKE SURE to change it in Particle.cpp (default constructor of a particle as well)
	hex.cleanFiles();

	for (int j = 0; j < NoOfEvolutions; j++)
	{

		// Refer to http://www.cplusplus.com/reference/ctime/time/ for below- done to get a different seed at every run
		time_t timer;
		time(&timer);  /* get current time; same as: timer = time(NULL); the number of seconds elapsed since 00:00 hours, Jan 1, 1970 UTC (i.e., a unix timestamp).   */
		struct tm y2k;
		double seconds;
		y2k.tm_hour = 5;   y2k.tm_min = 30; y2k.tm_sec = 0;
		y2k.tm_year = 75;	y2k.tm_mon = 0; y2k.tm_mday = 1;
		seconds = difftime(timer, mktime(&y2k));	/* the number of seconds elapsed since given by y2k time; tm_year = x means 1900+x years  */
		seed = -1 * timer;
		seed2 = -1 * seconds;
		//printf("%lld \t %.f \t %s \n",mktime(&y2k),seconds,asctime(&y2k));
		
	//	hex.setHexagonalLattice();
	//	hex.setConstraintCircle();	// NO USE
	//	hex.setConstraintHexagon();
	
		//hex.getNeighbours_ConstraintHex();		// To check initaial lattice and neighbours
			
	//  For Cell Position input from MATLAB
		cout << "Data Input from MATLAB" << endl;
		hex.InitialGeometry();
		hex.EvolutionType();

		hex.setParticlePosition_Matlab();
		hex.setNeighbourConnection_Matlab();
		hex.calculateEquilibriumDistance_Matlab();
		//cout << "Main.cpp y "; cin >> y;

		hex.test_checkNeighbours_ConstraintHex();	// can be used for input from MATLAB as well
	//	cout << "Test_Check Neighbours; Input any integer "; cin >> y;

		cout << "Total No. of Iterations in an Evolution " << totalIterations << endl;
		cout << "Total No. of Evolutions " << NoOfEvolutions << endl;
		hex.printParameters();

		hex.setRandomPolarization(seed);
	   // hex.setBoundaryPolarization();				// To set Polarization of particles at boundary
		hex.setBoundaryPolarization_Matlab();
	//	hex.setParticleParameters(seed2, j);
	//	hex.printLattice();						// prints particle characteristics on terminal at t=0
	//	cout << "Lattice Size " << latticeLength << "x" << latticeBreadth << endl;
			
		cout << "Simulation No " << j + 1 << endl;
		
		for (int i = 0; i < totalIterations+1; i++)
		{
	
		//	to print all the characteristics of a particle at a particular iteration
			/*if (i == 0 || i== 1)			
			{
				cout << "in check evolution loop at " <<i<<" iteration"<< endl;
				particleNo = 13;
				Particle_checkEvolution = hex.getParticle(particleNo);
				Particle_checkEvolution.print_checkEvolution(particleNo, i, count);   //count is used to keep track of how many no. of times this loop has been called
				ofstream out;
				out.open("checkEvolution.txt", ios::app);
				out << endl;
				out.close();
				count++;
			}*/

		//	function to write file for quiver Plot
			if (i % QuiverplotInterval == 0 || i == 1)
			{
				hex.forQuiverplot(i);
				cout << "writing to file for quiver plot at t = " << i << endl;
			}
			
			if (i == 0)
			{
				cout << "To continue with the above desired settings, type any integer" << endl;
				cin >> falseInput;
			}
	   	//  function to write file for force histogram plot   NOTE: ALWAYS CHECK THE CONDITION AT WHICH DATA POINTS ARE TAKEN BEFORE SIMULATING
			//if ( ((i%forceHistInterval==0) && (i>=100)) || i==500)
			//{
			//	cout << "writing to file for force histogram plot at t = " << i << endl;
			//	while (RowNo <= latticeBreadth)
			//	{
			//		//hex.ForceHistogramPlot(RowNo);
			//		hex.newForceHistogram(RowNo, i);
			//		count_RowNo++;
			//		RowNo = (StartRowNo - 1) + (count_RowNo - 1)*RowNoInterval;
			//	}
			//	RowNo = StartRowNo;
			//	count_RowNo = 1;
			//}

		//	function to write to file for velocity contour plot
			/*hex.VelocityContourPlot(i);		// giving iteration No as input
			if (i % 2000 == 0){cout<<i<<endl;}*/
	

		//	hex.Evolve_hex();  //for old evolution type
			
		
		//	hex.VelOrderParameter(j);
			// For Matlab Evolution
			//if (i % 100 == 0){ cout << "timestep "<< i << endl; }
			hex.Evolve_ConstraintHex_Final();   // when input is taken from Matlab
		
			//	hex.test_checkNeighbours_ConstraintHex();
			//cin >> y;

	    	//hex.appendFiles(NoOfEvolutions);
		}
	}
	cout << "Evolution complete. Type an integer to continue" << endl;
	cin >> y;
	return 0;
}
