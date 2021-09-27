class Position{
	public:
		double a;
		Position();		// needed to for neighbours array
		double rx;
		double ry;
};

class Velocity{
	public:
		Velocity();		// to set initial velocity to zero
		double vx;
		double vy;
		double modVel;
};

class Polarization{
	public:
		double px;
		double py;
		double modPol;
		Polarization();       
};

class Force
{
public:
	double fx;
	double fy;
	Force();
};

class Neighbours
{
public:
	Position pos;
	double EquilibriumDistance;
	//Neighbours();
};

class Noise
{
public:
	Noise();
	double nx;
	double ny;
};

class Particle
{
	Position pos;
	Velocity vel;
	Polarization pol;
	Force force;
	Force force_pol[2];	//force acting along ([0]  element) and perpendicular ([1] element) to polarization in x-y cordinates
	Velocity velocity_pol[2]; //velocity acting along ([0]th  element) and perpendicular ([1] element) to polarization in x-y cordinates
	Noise noise[2];  //1st for Velocity and second for Position

public:
	double h;	//delta-t (time interval) between 2 updates in euler method
	double mu1;
	double mu2;
	double zeta;
	double v0;
	double v0_Boundary;
	double k;
	double k_extension;
	double k_compression;
	double noisePol;  //Random noise for PolarizationEvolution
	double D[2];   //Noise strengths- [1] for velocity and [2] for Polarizaition
	double NoisePreFactor[2];   //Noise Prefactor- [1] for velocity and [2] for Polarizaition

	int count[2];	// to keep track of whether threshold for motility is crossed
	int boundary;	//to check if the paricle is at Boundary
	double  thresholdVel;
	int * NeighbourNo_Matlab;
	Neighbours * NeighbourArray;
	double ConfinementForce;

//	double lambda;

	Particle();
//	~Particle();

	void SetParameters();
	void GetParameters(double * Parameters);
	//double Getv0_Boundary();

	void SetPosition(Position p);
	void SetPosition(double x, double y);
	Position GetPosition();
	Position positionEvolution(int evolutionType);
	Position positionEvolution_Constraint(int evolutionType, double * LatticeParameters, int flag_Matlab, double * MatlabParameters_Circle);

	void SetVelocity(Velocity v);
	Velocity GetVelocity();

	void SetPolarization(double x, double y);
	Polarization GetPolarization();
	Polarization polarizationEvolution(int evolutionType);


	void calculateForce(Position * neighbours);
	void SetForce(double x, double y);
	Force GetForce();

	void  calculateForce_pol();
	void SetForce_pol(Force * Force_pol);
	Force * GetForce_pol();
		
	Velocity * calculateVelocity_pol();
	Velocity * GetVelocity_pol();

	void setNeighbours();
	void print_checkEvolution(int particleNo, int time_iteration, int count);

	void setNeighbourNo_Matlab(int NeighbourNo, int j);
	int getNeighbourNo_Matlab(int j);

	void SetNoise(Noise * noise);
	
	double DotProduct(double x1, double y1, double x2, double y2);
	double Modulus(double x, double y);

	void writeFile();

};
