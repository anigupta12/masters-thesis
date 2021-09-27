#include "Particle.h"

class Lattice{
	int L;
	int B;
	double a;
	int N; //Total no.of particles
	Particle * ParticleNo;
	//Particle * BoundaryParticle;
	//Particle ** ParticleArray;
	int evolutionType;		// 0 for simple evolution; 1 for evolution with motility of inside cells to be zero initially
	int flag_Matlab;       // To check if input is taken from Matlab
	int geometryType;		// To know the initial geometry- circular or rectangular mesh
	double MatlabParameters_circle[4];	//x0,y0 (Centre of Circle), Boundary/Confinement Radius, Hard Boundary (in order); Written in setParticlePosition_Matlab() function;
	int MaxNeighbours_Matlab;	//Maximum Possible Neighbours in Matlab		

public:
	Lattice(int l, int b, double t);
	void getLatticeParameters(double * LatticeParameters);
	void setHexagonalLattice();
	void setRandomPolarization(long seed);
	void setParticleParameters(long seed, int simulationNo);
	void setBoundaryPolarization();
	void setNoise();
	void getNeighbours_hex(int ParticleNo,Position * neighbours);
	void printLattice();						// to print initial lattice on terminal
	void printParameters();						// to print parameter values on terminal
	Particle getParticle(int ParticleNo);		// to get all the information of particle at time t without having to open the .txt file
	void EvolutionType();						// To get input of evolution type- 1 for normal and 2 for inner cells motility set to zero initially
	void Evolve_hex();
	double ran2(long *idum);	//random variable from Uniform distribution
	double gasdev(long *idum);  //random variable from Gaussian distribution
	void cleanFiles();
	void appendFiles(int EvolutionNo);
	void forQuiverplot(int timestep);
	void ForceHistogramPlot(int rowNo);
	void VelocityContourPlot(int iterationNo);
//  test functions
	void test_checkNeighbours(int ParticleNo, Position * neighbours);
	void test_calculateForce(int ParticleNo, Position * neighbours);
	
	void setConstraintCircle();
	void setConstraintHexagon();
	void getNeighbours_ConstraintHex();
	void calculateForce_ConstraintHex(int ParticleNo);
	//void Lattice::Evolve_ConstraintHex();
	void Evolve_ConstraintHex_Final();
	void test_checkNeighbours_ConstraintHex();

	void InitialGeometry();
	void setParticlePosition_Matlab();
	void setNeighbourConnection_Matlab();
	void getNeighbourConnection_Matlab();
	void calculateForce_Matlab(int i);
	void calculateEquilibriumDistance_Matlab();
	void setBoundaryPolarization_Matlab();

	void VelOrderParameter(int EvolutionNo);
	void newForceHistogram(int rowNo, int timestep);
	//void forQuiverplot1();
	//void forQuiverplot2();
};
