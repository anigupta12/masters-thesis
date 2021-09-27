//#define _CRT_SECURE_NO_DEPRECATE
#include "Lattice.h"
#include <iostream>
#include <cmath>
#include<fstream>
#include <sstream>
#include<stdio.h>
#include<time.h>

using namespace std;

Lattice::Lattice(int l, int b, double t)
{
	L = l;  //Total no. of particles in a row
	B = b;  //Total no. of particles in a column
	a = t;  //Distance between 2 points
	//N = L*B;
	//ParticleNo = new Particle[N];
	//BoundaryParticle = new Particle[2 * L + 2 * B];   //NO USE
	evolutionType = 0;
	flag_Matlab = 0;
}

void Lattice::getLatticeParameters(double * LatticeParameters)
{
	LatticeParameters[0] = L;
	LatticeParameters[1] = B;
	LatticeParameters[2] = a;
}

void Lattice::setHexagonalLattice()
{
	double i = 0;	// temp variable to increment column
	int j = 0;		// temp variable to increment row
	int n = 0;		// temp variable to count no. of particles
	double sin60 = sqrt(3) / 2;
	double cos60;
	cos60 = 0.5;
	ofstream out;
	out.open("lattice.dat");
	Position p;      //temporary variable to store updated position and then write to a file
	
	N = L*B;		//N is total no. of Particles
	ParticleNo = new Particle[N];

	while (j<B)		 // B is breadth of the lattice (no. of points in a column) and L is the lenght of the lattice (no. of points in a row)	
	{
		if (j % 2 == 0)
		{
			i = 0;
		//	printf("\n");
			while (i < L){
				ParticleNo[n].SetPosition(a*i, a*j*sin60);
				p = ParticleNo[n].GetPosition();
			//	printf("%.10e\t%.10e\n", p.rx, p.ry);
				out << p.rx << "\t" << p.ry << endl;
				i++;	//increment in column
				n++;	// to count no. of particles
			}
		}
		else
		{
			i = 0;
		//	printf("\n");
			while (i < L){
				ParticleNo[n].SetPosition(a*(i+cos60), a*j*sin60);
				p = ParticleNo[n].GetPosition();
			//	printf("%.10e\t%.10e\n", p.rx, p.ry);
				out << p.rx << "\t" << p.ry <<endl;
				i++;	//increment in column
				n++;		// to count no. of particles
			}
		}
		j++;	//increment in row
	}
	N = n;		// to count no. of particles
}

// void Lattice::setConstraintHexagon()     // Make sure B is odd
// {
// 	int n = 0;
// 	double x = 0, y = 0;
// 	double sin60 = sqrt(3) / 2;
// 	double cos60;
// 	cos60 = 0.5;

// 	N = L*B + (B - 1) / 2 * (B - 1) / 2;
// 	ParticleNo = new Particle[N];
	
// 	for (int i = 0; i < B; i++)
// 	{
// 		if (i <= (B - 1) / 2)
// 		{
// 			for (int j = 0; j < L+i; j++)
// 			{
// 				x = (B - 1) / 2 * a*cos60 + j*a - i*a*cos60;
// 				y = i*sin60*a;
// 				ParticleNo[n].SetPosition(x, y);
// 				n++;
// 				/*if (j==0 || j==L-1)
// 				{
// 				BoundaryParticle
// 				}*/
// 			//	cout << x<<" ";
// 			}
// 			//cout << endl;
// 		}
// 		else
// 		{
// 			for (int j = 0; j < (L + (B-1)/2 - (i-(B-1)/2)); j++)
// 			{
// 				x = j*a + (i - (B - 1) / 2)*a*cos60;
// 				y = i*sin60*a;
// 				ParticleNo[n].SetPosition(x, y);
// 				n++;
// 			//	cout << x<<" ";
// 			}
// 		//	cout << endl;
// 		}
// 	}
// 	N = n;
// 	cout << "Constraint Hexagonal Lattice. Total No. of Particles- "<<N << endl;;
// 	//cin >> n;
// }

void Lattice::setConstraintCircle()		// Make sure L is odd and B=2*L
{
	double sin60 = sqrt(3) / 2;
	double cos60;
	cos60 = 0.5;
	double x = 0, y = 0;
	Position temp;
	double CircleBoundary;
	//Particle tempa;
	for (int i = 0; i < N; i++)
	{
		temp=ParticleNo[i].GetPosition();
		// cout << temp.rx << "\t" << temp.ry << endl;
		// cin >> i;
		CircleBoundary = (temp.rx - L) * (temp.rx - L) + (temp.ry - L*sin60) *(temp.ry - L*sin60);
		if (CircleBoundary>L*L)
		{
			ParticleNo[i].SetPosition(-1*10*a, -1*10*a);
		}
	}
}

// For Input from MATLAB
void Lattice::InitialGeometry()
{
	cout << "Initial geometry of lattice: Enter 1 for Rectangular Mesh and 2 for Circular Confinement" << endl;
	cin >> geometryType;
}
void Lattice::setParticlePosition_Matlab()
{
	//flag_Matlab = 1;
	int i = 0;
	//int lSize;
	Position temp_pos;
	FILE * positionFile;
	FILE * pFile;
	positionFile = fopen("position_matlab.dat", "r");
	//fseek(position, 0, SEEK_END);
	//lSize = ftell(position);
    //rewind(position);
	//cout << lSize << endl;
	//cin >> lSize;
	//for (int i = 0; i < N; i++)
	pFile = fopen("position_matlab_check.dat", "w");
	ifstream fin;
	fin.open("position_matlab.dat");
	cout << "What are total no. of Particles ";
	cin >> N;
	ParticleNo = new Particle[N];
	while (i<N)//(fgetc(positionFile) != EOF)//(!feof(positionFile))
	{
		
		fin >> temp_pos.rx;
		fin >>temp_pos.ry;

		//fscanf(positionFile, "%e", &temp_pos.rx);
		//fscanf(positionFile, "%e", &temp_pos.ry);
		ParticleNo[i].SetPosition(temp_pos);
		
		// To check whether 
		fprintf(pFile, "%.10e\t%.10e\n", ParticleNo[i].GetPosition().rx, ParticleNo[i].GetPosition().ry);
		i++;
	}
	fclose(positionFile);
	fclose(pFile);
	fin.close();
	//cout << "Printing Total No. of Particles "<<N << endl;
	//cout << temp_pos.rx << " " << temp_pos.ry << endl;
	//cin >> i;
	if (geometryType == 2)
	{
		for (int k = 0; k < 4; k++)
		{
			if (k == 0)
			{
				cout << "Centre of Circle x0 and y0: " << endl;
			}
			if (k == 2)
			{
				cout << "Radius of Boundary/Inner Confinement and Hard Boundary/Outer confinement : " << endl;
			}
			cin >> MatlabParameters_circle[k];
		}
	}
	cout << "Maximum Possible Neighbours possible " << endl;
	cin >> MaxNeighbours_Matlab;
}
void Lattice::setNeighbourConnection_Matlab()
{
		flag_Matlab = 1;
		int i = 0;
		int random = 0;
		int MaxNeighbour = MaxNeighbours_Matlab;
		int NeighbourNo;
		//NeighbourNo = new int [MaxNeighbour];
		FILE * positionFile;
		FILE * pFile;
		positionFile = fopen("neighbour_matlab.dat", "r");
		pFile = fopen("neighbour_matlab_check.dat", "w");
		ifstream fin;
		fin.open("neighbour_matlab.dat");
		while (i<N)//(fgetc(positionFile) != EOF)//(!feof(positionFile))
		{
			for (int j = 0; j < MaxNeighbour; j++)
			{
				fin >> NeighbourNo;
				ParticleNo[i].setNeighbourNo_Matlab(NeighbourNo, j);
				fprintf(pFile, "%i \t", ParticleNo[i].getNeighbourNo_Matlab(j)); //To check whether Neighbour Connections are fine

				if (NeighbourNo != 0)
				{
					ParticleNo[i].NeighbourArray[j].pos = ParticleNo[ParticleNo[i].getNeighbourNo_Matlab(j) - 1].GetPosition();
					//if (j == 6 && i == 18){ cout << ParticleNo[i].NeighbourArray[j].pos.rx << " " << ParticleNo[i].NeighbourArray[j].pos.ry << endl; cout << i << " " << j << " " << NeighbourNo[j] - 1; cin>>random; }
				}
				//else  // To check non existing neighbours are ignored or not 
				//{
				//	ParticleNo[i].NeighbourArray[j].pos; 
				//	cout << ParticleNo[i].NeighbourArray[j].pos.rx << " " << ParticleNo[i].NeighbourArray[j].pos.ry << endl; cout << i << " " << j << " " << NeighbourNo[j] - 1;
				// cin >> random; }
			}
			fprintf(pFile, "\n");
			i++;
		}
		
		fclose(positionFile);
		fclose(pFile);
		fin.close();
}
void Lattice::getNeighbourConnection_Matlab()
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < MaxNeighbours_Matlab; j++)
		{
			if (ParticleNo[i].getNeighbourNo_Matlab(j) != 0)
			{
				ParticleNo[i].NeighbourArray[j].pos = ParticleNo[ParticleNo[i].getNeighbourNo_Matlab(j) - 1].GetPosition();
				//if (j == 6 && i == 18){ cout << ParticleNo[i].NeighbourArray[j].pos.rx << " " << ParticleNo[i].NeighbourArray[j].pos.ry << endl; cout << i << " " << j << " " << NeighbourNo[j] - 1; cin>>random; }
			}
			//else  // To check non existing neighbours are ignored or not 
			//{
			//	ParticleNo[i].NeighbourArray[j].pos; 
			//	cout << ParticleNo[i].NeighbourArray[j].pos.rx << " " << ParticleNo[i].NeighbourArray[j].pos.ry << endl; cout << i << " " << j << " " << NeighbourNo[j] - 1;
			// cin >> random; }
		}
	}
}
void Lattice::setBoundaryPolarization_Matlab()
{
	//int i; //temp variable for i-th particle
	if (geometryType == 1)
	{
		int temp_particleNo;
		ifstream fin1;
		ifstream fin2;
		fin1.open("topBoundary_matlab.dat");
		fin2.open("bottomBoundary_matlab.dat");
		while (!fin1.eof())
		{
			fin1 >> temp_particleNo;
			ParticleNo[temp_particleNo - 1].boundary = 1;
			ParticleNo[temp_particleNo - 1].SetPolarization(0, 1);
			ParticleNo[temp_particleNo - 1].v0 = ParticleNo[temp_particleNo - 1].v0_Boundary;
			ParticleNo[temp_particleNo - 1].count[0] = 1;
			ParticleNo[temp_particleNo - 1].count[1] = 2;

		}
		while (!fin2.eof())
		{
			fin2 >> temp_particleNo;
			ParticleNo[temp_particleNo - 1].boundary = 1;
			ParticleNo[temp_particleNo - 1].SetPolarization(0, -1);
			ParticleNo[temp_particleNo - 1].v0 = ParticleNo[temp_particleNo - 1].v0_Boundary;
			ParticleNo[temp_particleNo - 1].count[0] = 1;
			ParticleNo[temp_particleNo - 1].count[1] = 2;
		}
		fin1.close();
		fin2.close();
	}
}

void Lattice::setRandomPolarization(long seed)
{
	int i; //temp variable for ith particle
	double theta;
	for ( i = 0; i < N; i++)
	{
		theta = 8 * atan(1.0)*(ran2(&seed) - 0.5);
		ParticleNo[i].SetPolarization(cos(theta), sin(theta));
	}
}

void Lattice::setBoundaryPolarization()
{
	int i; //temp variable for i-th particle
	for (i = 0; i < L; i++){
		ParticleNo[i].SetPolarization(0, -1);
		ParticleNo[i].v0 = 1.5;
		ParticleNo[i].count[0] = 1;
		ParticleNo[i].count[1] = 2;
		}
	for (i = L*(B - 1); i <N ; i++)
	{
		ParticleNo[i].SetPolarization(0, 1);
		ParticleNo[i].v0 = 1.5;
		ParticleNo[i].count[0] = 1;
		ParticleNo[i].count[1] = 2;
	}
}

void Lattice::EvolutionType()
{
	cout << "To set the motility of inner cells to be zero (v0=0) initially, enter 1 else 0" << endl;
	cin >> evolutionType;
}

void Lattice::setParticleParameters(long seed, int simulationNo)
{
	cout << "Individual Particle parameters are being set" << endl;
	int count[2];
	count[2] = 0;

//	setting inhomogeneous spring stiffness
	double lowLimit, highLimit;
	if (simulationNo == 0)
	{
		cout << "For inhomogeneous spring stiffness enter 1 else 0" << endl;
		cin >> count[0];
	}
	if (count[0] == 1){
		cout << "Enter the range for spring stiffness (lower limit and higher limit)" << endl;
		cin >> lowLimit >> highLimit;
		cout << "Inhomogeneous spring stiffness uniformly distributed from [" << lowLimit << " " << highLimit << "]" << endl<<endl;
		for (int i = 0; i < N; i++)
		{
			ParticleNo[i].k = (highLimit-lowLimit) * ran2(&seed) + lowLimit;
		}
	}

//	setting motility (and polarization) of inner cells to be zero (v0=0)
	if (simulationNo==0)
	{
		cout << "For setting motility (and polarization) of inner cells to be zero (v0=0) enter 1 else 0" << endl;
		cin >> count[1];
	}
	if (count[1] == 1){
		cout << "Motility (v0) and Polarization of inner cells set to 0 " << endl;
		evolutionType = 1;
		for (int i = 0; i < L; i++){
			ParticleNo[i].v0 = 0;
			ParticleNo[i].SetPolarization(0, 0);
		}
		for (int i = L*(B - 1); i < N; i++)
		{
			ParticleNo[i].v0 = 0;
			ParticleNo[i].SetPolarization(0, 0);
		}
	}
}

void Lattice::getNeighbours_hex(int n, Position * neighbours)	// n is the particle no. for which neighbours function is called
{
	int temp_i = 0;
	
	// if the particle is on extreme left
	if (n%L==0)
	{
		if (n==0)		// for particle at bottom extreme left
		{
			neighbours[0] = ParticleNo[n+1].GetPosition();
			neighbours[1] = ParticleNo[n+L].GetPosition();
			neighbours[2] = ParticleNo[n + L - 1].GetPosition();    //cylindrical boundary condition
			neighbours[2].rx = neighbours[2].rx - L*a;
			neighbours[3] = ParticleNo[n + (2 * L) - 1].GetPosition(); //cylindrical boundary condition
			neighbours[3].rx = neighbours[3].rx - L*a;
		}
		else if (n==L*(B-1))			//for particle at top extreme left
		{
			temp_i = L*(B - 1);
			if (B%2 == 0)
			{
				neighbours[0] = ParticleNo[temp_i - (L - 1)].GetPosition();
				neighbours[1] = ParticleNo[temp_i - L].GetPosition();
				neighbours[2] = ParticleNo[temp_i + 1].GetPosition();
				neighbours[3] = ParticleNo[temp_i + (L - 1)].GetPosition();  //cylindrical boundary condition
				neighbours[3].rx = neighbours[3].rx - L*a;
			}
			else
			{
				neighbours[0] = ParticleNo[temp_i - L].GetPosition();
				neighbours[1] = ParticleNo[temp_i + 1].GetPosition();
				neighbours[2] = ParticleNo[temp_i + (L - 1)].GetPosition();  //cylindrical boundary condition
				neighbours[2].rx = neighbours[2].rx - L*a;
				neighbours[3] = ParticleNo[temp_i - 1].GetPosition();		//cylindrical boundary condition
				neighbours[3].rx = neighbours[3].rx - L*a;
			}
		}
		else if (n%(2 * L) == 0)     // for rest of the particles on left side: there will be 2 conditions
		{
			neighbours[0] = ParticleNo[n+1].GetPosition();
			neighbours[1] = ParticleNo[n-L].GetPosition();
			neighbours[2] = ParticleNo[n+L].GetPosition();
			neighbours[3] = ParticleNo[n + L - 1].GetPosition(); //cylindrical boundary condition
			neighbours[3].rx = neighbours[3].rx - L*a;
			neighbours[4] = ParticleNo[n + (2*L) - 1].GetPosition(); //cylindrical boundary condition
			neighbours[4].rx = neighbours[4].rx - L*a;
			neighbours[5] = ParticleNo[n - 1].GetPosition(); //cylindrical boundary condition
			neighbours[5].rx = neighbours[5].rx - L*a;
		}
		else
		{
			neighbours[0] = ParticleNo[n + 1].GetPosition();
			neighbours[1] = ParticleNo[n - L].GetPosition();
			neighbours[2] = ParticleNo[n + L].GetPosition();
			neighbours[3] = ParticleNo[n + L - 1].GetPosition(); //cylindrical boundary condition
			neighbours[3].rx = neighbours[3].rx - L*a;
			neighbours[4] = ParticleNo[n + L + 1].GetPosition();
			neighbours[5] = ParticleNo[n - L + 1].GetPosition();
		}
	}

	// if the particle is on extreme right
	else if ((n+1)%L==0)
	{
		if (n==(L-1))		// for particle at bottom right
		{
			neighbours[0] = ParticleNo[n - L + 1].GetPosition();	//cylindrical boundary condition
			neighbours[0].rx = neighbours[0].rx + L*a;
			neighbours[1] = ParticleNo[n + L].GetPosition();
			neighbours[2] = ParticleNo[n + L - 1].GetPosition();
			neighbours[3] = ParticleNo[n - 1].GetPosition();
		}
		else if (n==(L*B-1))	// for top right
		{
			temp_i = L*B - 1;
			if (B%2==0)
			{
				neighbours[0] = ParticleNo[temp_i - L + 1].GetPosition();	//cylindrical boundary condition
				neighbours[0].rx = neighbours[0].rx + L*a;
				neighbours[1] = ParticleNo[temp_i - L].GetPosition();
				neighbours[2] = ParticleNo[temp_i - 1].GetPosition();
				neighbours[3] = ParticleNo[temp_i - (2 * L) + 1].GetPosition();	//cylindrical boundary condition
				neighbours[3].rx = neighbours[3].rx + L*a;
			}
			else
			{
				neighbours[0] = ParticleNo[temp_i - L + 1].GetPosition();	//cylindrical boundary condition
				neighbours[0].rx = neighbours[0].rx + L*a;
				neighbours[1] = ParticleNo[temp_i - L].GetPosition();
				neighbours[2] = ParticleNo[temp_i - 1].GetPosition();
				neighbours[3] = ParticleNo[temp_i - L - 1].GetPosition();
			}
		}
		else if ((n+1)%(2 * L) == 0)  // for rest of the particles on right side: there will be 2 conditions
		{
			neighbours[0] = ParticleNo[n - 1].GetPosition();
			neighbours[1] = ParticleNo[n - L].GetPosition();
			neighbours[2] = ParticleNo[n + L].GetPosition();
			neighbours[3] = ParticleNo[n - L + 1].GetPosition(); //cylindrical boundary condition
			neighbours[3].rx = neighbours[3].rx + L*a;
			neighbours[4] = ParticleNo[n + 1].GetPosition(); //cylindrical boundary condition
			neighbours[4].rx = neighbours[4].rx + L*a;
			neighbours[5] = ParticleNo[n - (2 * L) + 1].GetPosition(); //cylindrical boundary condition
			neighbours[5].rx = neighbours[5].rx + L*a;
		}
		else
		{
			neighbours[0] = ParticleNo[n - 1].GetPosition();
			neighbours[1] = ParticleNo[n - L].GetPosition();
			neighbours[2] = ParticleNo[n + L].GetPosition();
			neighbours[3] = ParticleNo[n - L + 1].GetPosition(); //cylindrical boundary condition
			neighbours[3].rx = neighbours[3].rx + L*a;
			neighbours[4] = ParticleNo[n + L - 1].GetPosition();
			neighbours[5] = ParticleNo[n - L - 1].GetPosition();
		}
	}
	else if (n<L)		// for the particles in lowest line 
	{
		neighbours[0] = ParticleNo[n - 1].GetPosition();
		neighbours[1] = ParticleNo[n + 1].GetPosition();
		neighbours[2] = ParticleNo[n + L].GetPosition();
		neighbours[3] = ParticleNo[n + L - 1].GetPosition();
	}
	else if (n > L*(B - 1))		// for the particles on topmost line
	{
		if (B % 2 == 0)
		{
			neighbours[0] = ParticleNo[n - 1].GetPosition();
			neighbours[1] = ParticleNo[n + 1].GetPosition();
			neighbours[2] = ParticleNo[n - L].GetPosition();
			neighbours[3] = ParticleNo[n - L + 1].GetPosition();
		}
		else
		{
			neighbours[0] = ParticleNo[n - 1].GetPosition();
			neighbours[1] = ParticleNo[n + 1].GetPosition();
			neighbours[2] = ParticleNo[n - L].GetPosition();
			neighbours[3] = ParticleNo[n - L - 1].GetPosition();
		}
	}
	else						// for rest of the particles in centre
	{
		temp_i = n/L;
		if (temp_i % 2 == 0)
		{
			neighbours[0] = ParticleNo[n - 1].GetPosition();
			neighbours[1] = ParticleNo[n + 1].GetPosition();
			neighbours[2] = ParticleNo[n + L].GetPosition();
			neighbours[3] = ParticleNo[n - L].GetPosition(); 
			neighbours[4] = ParticleNo[n + L - 1].GetPosition();
			neighbours[5] = ParticleNo[n - L - 1].GetPosition();
		}
		else
		{
			neighbours[0] = ParticleNo[n - 1].GetPosition();
			neighbours[1] = ParticleNo[n + 1].GetPosition();
			neighbours[2] = ParticleNo[n + L].GetPosition();
			neighbours[3] = ParticleNo[n - L].GetPosition();
			neighbours[4] = ParticleNo[n + L + 1].GetPosition();
			neighbours[5] = ParticleNo[n - L + 1].GetPosition();
		}
	}
}

void Lattice::getNeighbours_ConstraintHex()  //Make sure that B=2*L-1
{
	int temp_i = 0;
	int n = 0;
	int temp_L = L;
	//double x = 0, y = 0;
	//double sin60 = sqrt(3) / 2;
	//double cos60;
	//cos60 = 0.5;
	for (int i = 0; i < B; i++)
	{
		if (i < (B - 1) / 2)
		{
			for (int j = 0; j < L+i; j++)
			{
				if (i == 0)		// temp_L=L for this case
				{
					if (j == 0)
					{
						ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
						ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n + L].GetPosition();
						ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n + L + 1].GetPosition();
						//ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n + L - 1].GetPosition();  //cylindrical boundary condition
						//ParticleNo[n].NeighbourArray[3].pos.rx = ParticleNo[n].NeighbourArray[3].pos.rx - L*a;
						n++;
					}
					else if (j == L - 1)
					{
						ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + L].GetPosition();
						ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n + L + 1].GetPosition();
						ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n - 1].GetPosition();
						//ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n - L + 1].GetPosition();  //cylindrical boundary condition
						//ParticleNo[n].NeighbourArray[3].pos.rx = ParticleNo[n].NeighbourArray[3].pos.rx + L*a;
						n++;
						temp_L = temp_L + 1;
			//			cout << "1 "<<temp_L << endl;
					}
					else
					{
						ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
						ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n - 1].GetPosition();
						ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n + L].GetPosition();
						ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n + L + 1].GetPosition();
						n++;
					}
				}
				else
				{
					if (j == 0)
					{
						ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
						ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n + temp_L].GetPosition();
						ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n + temp_L + 1].GetPosition();
						ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n - temp_L + 1].GetPosition();
						//ParticleNo[n].NeighbourArray[4].pos = ParticleNo[n + temp_L - 1].GetPosition();  //cylindrical boundary condition
						//ParticleNo[n].NeighbourArray[4].pos.rx = ParticleNo[n].NeighbourArray[4].pos.rx - L*a;
						n++;
					}
					else if (j == L+i - 1)
					{
						ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + temp_L].GetPosition();
						ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n + temp_L + 1].GetPosition();
						ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n - 1].GetPosition();
						ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n - temp_L].GetPosition();
						//ParticleNo[n].NeighbourArray[4].pos = ParticleNo[n - temp_L + 1].GetPosition();  //cylindrical boundary condition
						//ParticleNo[n].NeighbourArray[4].pos.rx = ParticleNo[n].NeighbourArray[4].pos.rx + L*a;
						n++;
						temp_L = temp_L + 1;
			//			cout << "2 "<<temp_L << endl;
					}
					else
					{
						ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
						ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n - 1].GetPosition();
						ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n + temp_L].GetPosition();
						ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n + temp_L + 1].GetPosition();
						ParticleNo[n].NeighbourArray[4].pos = ParticleNo[n - temp_L].GetPosition();
						ParticleNo[n].NeighbourArray[5].pos = ParticleNo[n - temp_L + 1].GetPosition();
						n++;
					}
				}
			}
		}
		if (i==(B-1)/2)
		{
			for (int j = 0; j < L+i; j++)
			{
				if (j == 0)
				{
					ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
					ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n + temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n - temp_L + 1].GetPosition();
					// No cylindrical Boundary condition
					n++;
				}
				else if (j == L+i - 1)
				{
					ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + temp_L - 1].GetPosition();
					ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n - 1].GetPosition();
					ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n - temp_L].GetPosition();
					// No cylindrical Boundary condition
					n++;
					temp_L = temp_L - 1;
			//		cout << "3 "<<temp_L << endl;
				}
				else
				{
					ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
					ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n - 1].GetPosition();
					ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n + temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n + temp_L - 1].GetPosition();
					ParticleNo[n].NeighbourArray[4].pos = ParticleNo[n - temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[5].pos = ParticleNo[n - temp_L + 1].GetPosition();
					n++;
				}
			}
		}
		if ( (i>(B-1)/2) && (i!=(B-1)) )
		{
			for (int j = 0; j < (L + (B - 1) / 2 - (i - (B - 1) / 2)); j++)
			{
				if (j == 0)
				{
					ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
					ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n + temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n - temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n - temp_L - 1].GetPosition();
					// No cylindrical Boundary condition
					n++;
				}
				else if (j == (L + (B - 1) / 2 - (i - (B - 1) / 2)) - 1)
				{
					ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + temp_L - 1].GetPosition();
					ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n - 1].GetPosition();
					ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n - temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n - temp_L - 1].GetPosition();
					// No cylindrical Boundary condition
					n++;
					temp_L = temp_L - 1;
			//		cout << "4 "<<temp_L<< endl;
				}
				else
				{
					ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
					ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n - 1].GetPosition();
					ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n + temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n + temp_L - 1].GetPosition();
					ParticleNo[n].NeighbourArray[4].pos = ParticleNo[n - temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[5].pos = ParticleNo[n - temp_L - 1].GetPosition();
					n++;
				}
			}
		}
		if (i==(B-1))
		{
			for (int j = 0; j < L; j++)
			{
				if (j == 0)
				{
					ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
					ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n - temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n - temp_L - 1].GetPosition();
					// No cylindrical Boundary condition
					n++;
				}
				else if (j == L - 1)
				{
					ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n - temp_L - 1].GetPosition();
					ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n - 1].GetPosition();
					ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n - temp_L].GetPosition();
					// No cylindrical Boundary condition
					n++;
			//		cout << "5 "<<temp_L << endl;
				}
				else
				{
					ParticleNo[n].NeighbourArray[0].pos = ParticleNo[n + 1].GetPosition();
					ParticleNo[n].NeighbourArray[1].pos = ParticleNo[n - 1].GetPosition();
					ParticleNo[n].NeighbourArray[2].pos = ParticleNo[n - temp_L].GetPosition();
					ParticleNo[n].NeighbourArray[3].pos = ParticleNo[n - temp_L - 1].GetPosition();
					n++;
				}
			}
		}
	}
	//cout << "In get Neighbours loop, N is "<<n << endl;
	/*cout << temp_L << endl;
	cin >> temp_i;*/
}

void Particle::calculateForce(Position * neighbours)		//To calculate force on a particle
{
	int i = 0;
	Position diff;		//difference between the particle and its neighbours
	double diff_mod;
	double temp_fx=0, temp_fy=0;	  //temporary variables to calculate force
	for (i = 0; i<6 ; i++)		
	{
		if (neighbours[i].rx != (-1 * 10 * pos.a))  // so that loop runs for specific no. of times according to no. of neighbours
		{
			diff.rx = neighbours[i].rx - pos.rx;
			diff.ry = neighbours[i].ry - pos.ry;
			diff_mod = Modulus(diff.rx, diff.ry);
			temp_fx = temp_fx + k*(diff_mod - pos.a)*(diff.rx / diff_mod);
			temp_fy = temp_fy + k*(diff_mod - pos.a)*(diff.ry / diff_mod);
		}
	}
	SetForce(temp_fx, temp_fy);
}

void Particle::calculateForce_pol()
{
	double temp_fx_pol = 0, temp_fy_pol = 0;		//temp_fx is modulus of force parallel to Polarization and temp_fy is perpendicular
	temp_fx_pol = DotProduct(pol.px, pol.py, force.fx, force.fy);	//mod of component of force along polarization & mod(Pol)=1
	force_pol[0].fx = temp_fx_pol*pol.px;		//component of force along polarization- in x-direction
	force_pol[0].fy = temp_fx_pol*pol.py;		//component of force along polarization- in y-direction
	force_pol[1].fx = force.fx - temp_fx_pol*pol.px;
	force_pol[1].fy = force.fy - temp_fx_pol*pol.py;
}

Position Particle::positionEvolution(int evolutionType)
{
	double force_alongPol;		// modulus of component of force parallel to polarization
	//double modVel;				// magnitude of velocity
	double temp_v0;
	double temp_mu2;
	force_alongPol = DotProduct(pol.px, pol.py, force.fx, force.fy);
	vel.modVel = Modulus(vel.vx, vel.vy);
	if (evolutionType==0)
	{
		vel.vx = (mu1*force_alongPol + v0)*pol.px + mu2*force_pol[1].fx + D[0] * noise[0].nx;
		vel.vy = (mu1*force_alongPol + v0)*pol.py + mu2*force_pol[1].fy + D[0] * noise[0].ny;
		pos.rx = pos.rx + h*vel.vx + NoisePreFactor[0] * noise[1].nx;
		pos.ry = pos.ry + h*vel.vy + NoisePreFactor[0] * noise[1].ny;
	}
	if (evolutionType==1)
	{
		if (vel.modVel>=thresholdVel)
		{
			count[0] = 1;
			count[1] = count[1]+1;
		}
		if (count[0]==0)
		{	
			temp_v0 = 0;
			temp_mu2 = mu1;
		}
		else
		{	
			temp_v0 = v0;
			temp_mu2 = mu2;
		}
		vel.vx = (mu1*force_alongPol + temp_v0)*pol.px + temp_mu2*force_pol[1].fx + D[0] * noise[0].nx;
		vel.vy = (mu1*force_alongPol + temp_v0)*pol.py + temp_mu2*force_pol[1].fy + D[0] * noise[0].ny;
		pos.rx = pos.rx + h*vel.vx + NoisePreFactor[0] * noise[1].nx;
		pos.ry = pos.ry + h*vel.vy + NoisePreFactor[0] * noise[1].ny;
	}
	return pos;
}

Polarization Particle::polarizationEvolution(int evolutionType)
{
	if (evolutionType==0)
	{
		pol.px = pol.px + h*(zeta*force_pol[1].fx) + NoisePreFactor[1] * noisePol * (-1) * pol.py;
		pol.py = pol.py + h*(zeta*force_pol[1].fy) + NoisePreFactor[1] * noisePol * 1 * pol.px;

	}
	if (evolutionType==1)
	{
		if (count[1] == 1)
		{
			pol.px = vel.vx / vel.modVel;
			pol.py = vel.vy / vel.modVel;
		}
		else
		{
			pol.px = pol.px + h*(zeta*force_pol[1].fx) + NoisePreFactor[1] * noisePol * (-1) * pol.py;
			pol.py = pol.py + h*(zeta*force_pol[1].fy) + NoisePreFactor[1] * noisePol * 1 * pol.px;
		}
	}
	pol.modPol = Modulus(pol.px, pol.py);
	pol.px = pol.px / pol.modPol;
	pol.py = pol.py / pol.modPol;
	return pol;
}

void Lattice::Evolve_hex()
{
	//Position check_InBoundary;
	for (int i = 0; i < N; i++)
	{
		Position neighbours[6];
		//check_InBoundary.rx = ParticleNo[i].GetPosition().rx;
		/*if (check_InBoundary.rx != -1 * 10 * a)
		{*/
		getNeighbours_hex(i, neighbours);
		ParticleNo[i].calculateForce(neighbours);
		ParticleNo[i].calculateForce_pol();
		//		ParticleNo[i].writeFile();
		//		test_checkNeighbours(i, neighbours);	// to check neighbours
		//		test_calculateForce(i, neighbours);		//	to check calculateForce f:n
		}
//	}
	for (int i = 0; i < N; i++)
	{
	/*	check_InBoundary.rx = ParticleNo[i].GetPosition().rx;
		if (check_InBoundary.rx != -1 * 10 * a)
		{*/
		ParticleNo[i].positionEvolution(evolutionType);
		ParticleNo[i].polarizationEvolution(evolutionType);
	  //}
	}
}

void Lattice::calculateForce_ConstraintHex(int n)		//To calculate force on a particle
{
	int i = 0;
	int TotalNeighbour;
	Position current;	// current position of the particle	
	Position diff;		//difference between the particle and its neighbours
	double diff_mod;
	double temp_fx = 0, temp_fy = 0;	  //temporary variables to calculate force
	current = ParticleNo[n].GetPosition();
	if (flag_Matlab==0)
	{
		TotalNeighbour = 6;
	}
	else if (flag_Matlab==1)
	{
		TotalNeighbour = MaxNeighbours_Matlab;
		//cout << "Total Neighbours are " << TotalNeighbour << endl;
	}
	for (i = 0; i<TotalNeighbour; i++)
	{
		if (ParticleNo[n].NeighbourArray[i].pos.rx != (-1 * 10 * a))  // so that loop runs for specific no. of times according to no. of neighbours
		{
			diff.rx = ParticleNo[n].NeighbourArray[i].pos.rx - current.rx;
			diff.ry = ParticleNo[n].NeighbourArray[i].pos.ry - current.ry;
			diff_mod = ParticleNo[n].Modulus(diff.rx, diff.ry);
			temp_fx = temp_fx + ParticleNo[n].k*(diff_mod - a)*(diff.rx / diff_mod);
			temp_fy = temp_fy + ParticleNo[n].k*(diff_mod - a)*(diff.ry / diff_mod);
		}
	}
	ParticleNo[n].SetForce(temp_fx, temp_fy);
}

void Lattice::calculateEquilibriumDistance_Matlab()
{
	int i = 0;
	int TotalNeighbour;
	Position current;	// current position of the particle	
	Position diff;		//difference between the particle and its neighbours
//	double diff_mod;
	FILE * pfile;
	pfile = fopen("EquilibriumDistance_matlab_check.dat", "w");
	TotalNeighbour = MaxNeighbours_Matlab;
	//cout << "Total Neighbours are " << TotalNeighbour << endl;
	int y;
	for (int n = 0; n < N; n++)
	{
		current = ParticleNo[n].GetPosition();
		for (i = 0; i < TotalNeighbour; i++)
		{
			if (ParticleNo[n].NeighbourArray[i].pos.rx != (-1 * 10 * a))  // so that loop runs for specific no. of times according to no. of neighbours
			{
				diff.rx = ParticleNo[n].NeighbourArray[i].pos.rx - current.rx;
				diff.ry = ParticleNo[n].NeighbourArray[i].pos.ry - current.ry;
				ParticleNo[n].NeighbourArray[i].EquilibriumDistance = ParticleNo[n].Modulus(diff.rx, diff.ry);
				fprintf(pfile, "%.10e\t", ParticleNo[n].NeighbourArray[i].EquilibriumDistance);
				/*if (i==2 || i==3)
				{*/
					//cout << ParticleNo[n].NeighbourArray[i].pos.rx << " " << current.rx << " " << diff.rx << endl;
					//cin >> y;
				//}
			}
		}
		fprintf(pfile, "\n");
	}
	fclose(pfile);
	
	cout << "Calculating Equillibrium Distance between two cells (natural spring length at t=0) " << endl;
}

void Lattice::calculateForce_Matlab(int n)
{
	int i = 0;
	int TotalNeighbour;
	Position current;	// current position of the particle	
	Position diff;		//difference between the particle and its neighbours
	double diff_mod;
	double springConstant;
	double x;           // spring compression/extension
	double EqDistance;	//Equilibrium Distance between particles at t=0
	double temp_fx = 0, temp_fy = 0;	  //temporary variables to calculate force
	current = ParticleNo[n].GetPosition();
	if (flag_Matlab == 0)
	{
		TotalNeighbour = 6;
	}
	else if (flag_Matlab == 1)
	{
		TotalNeighbour = MaxNeighbours_Matlab;
		//cout << "Total Neighbours are " << TotalNeighbour << endl;
	}
	for (i = 0; i<TotalNeighbour; i++)
	{
		if (ParticleNo[n].NeighbourArray[i].pos.rx != (-1 * 10 * a))  // so that loop runs for specific no. of times according to no. of neighbours
		{
			EqDistance = ParticleNo[n].NeighbourArray[i].EquilibriumDistance;
			diff.rx = ParticleNo[n].NeighbourArray[i].pos.rx - current.rx;
			diff.ry = ParticleNo[n].NeighbourArray[i].pos.ry - current.ry;
			diff_mod = ParticleNo[n].Modulus(diff.rx, diff.ry);
			x = diff_mod - EqDistance;
			if (x>0){ springConstant = ParticleNo[n].k_extension; }
			else { springConstant = ParticleNo[n].k_compression; }
			temp_fx = temp_fx + springConstant/EqDistance * x * (diff.rx / diff_mod);
			temp_fy = temp_fy + springConstant/EqDistance * x * (diff.ry / diff_mod);
		}
	}

	if (geometryType==2)		//Constraint Circle
	{
		double BoundaryRadius;			   // radius of confinement/ECM	
		//double radius_outer;			   // Hard Boundary
		double x0;
		double y0;
		double diff_x;
		double diff_y;
		double cos_theta;
		double sin_theta;
		double diff_modulus;
		double extension;
		//double flag_BoundaryCross = 0;	 // to check whether particle has crossed the boundary

		x0 = MatlabParameters_circle[0];
		y0 = MatlabParameters_circle[1];
		BoundaryRadius = MatlabParameters_circle[2];
		diff_x = (current.rx - x0);
		diff_y = (current.ry - y0);
		if (diff_x == 0 || diff_y == 0){ diff_modulus = 1; }
		else{ diff_modulus = ParticleNo[n].Modulus(diff_x, diff_y); }
		cos_theta = diff_x / diff_modulus;
		sin_theta = diff_y / diff_modulus;
		extension = diff_modulus - BoundaryRadius;
		if (extension>0)
		{
			temp_fx = temp_fx - ParticleNo[n].ConfinementForce*extension*cos_theta;
			temp_fy = temp_fy - ParticleNo[n].ConfinementForce*extension*sin_theta;
			//flag_BoundaryCross = 1.0;
			/*cout << "In confinement force condition" << endl;
			cin>>temp_cout;*/
		}
	}
	ParticleNo[n].SetForce(temp_fx, temp_fy);
}

void Lattice::setNoise()
{
	long seed1;
	long seed2;
	long seed3;
	long seed4;
	long seed5;

	// Refer to http://www.cplusplus.com/reference/ctime/time/ for below- done to get a different seed at every run
	time_t timer;
	time(&timer);  /* get current time; same as: timer = time(NULL); the number of seconds elapsed since 00:00 hours, Jan 1, 1970 UTC (i.e., a unix timestamp).   */
	struct tm y2k;
	double seconds;
	y2k.tm_hour = 6;   y2k.tm_min = 30; y2k.tm_sec = 0;
	y2k.tm_year = 75;	y2k.tm_mon = 0; y2k.tm_mday = 1;
	seconds = difftime(timer, mktime(&y2k));	/* the number of seconds elapsed since given by y2k time; tm_year = x means 1900+x years  */
	seed1 = -1 * seconds;
	y2k.tm_mday = 2;
	seconds = difftime(timer, mktime(&y2k));
	seed2 = -1 * seconds;
	y2k.tm_mday = 3;
	seconds = difftime(timer, mktime(&y2k));
	seed3 = -1 * seconds;
	y2k.tm_mday = 4;
	seconds = difftime(timer, mktime(&y2k));
	seed4 = -1 * seconds;
	y2k.tm_mday = 5;
	seconds = difftime(timer, mktime(&y2k));
	seed5 = -1 * seconds;

	for (int i = 0; i < N; i++)
	{
		Noise temp_noise[2];
		temp_noise[0].nx = gasdev(&seed1);	//noise for velocity (x-direction)
		temp_noise[0].ny = gasdev(&seed2);  //noise for velocity (y-direction)
		temp_noise[1].nx = gasdev(&seed3);  //noise for position (x-direction)
		temp_noise[1].ny = gasdev(&seed4);  //noise for position (y-direction) 
		ParticleNo[i].SetNoise(temp_noise);	
		ParticleNo[i].noisePol = gasdev(&seed5);
	}
}

void Lattice::Evolve_ConstraintHex_Final()
{
	if (flag_Matlab == 0)
	{
		getNeighbours_ConstraintHex(); 
	}
	else if (flag_Matlab == 1)
	{ 
		getNeighbourConnection_Matlab(); 
		//cout << "flag_Matlab=1 and Neighbour connection " << endl;
		//test_checkNeighbours_ConstraintHex();	// can be used for input from MATLAB as well
		//int y;
		//cout << "Test_Check Neighbours; Input any integer "; cin >> y;
	}
	for (int i = 0; i < N; i++)
	{
		if (flag_Matlab == 0)
		{ 
			calculateForce_ConstraintHex(i); 
		}
		else if (flag_Matlab==1)
		{
			calculateForce_Matlab(i); 
		}
		ParticleNo[i].calculateForce_pol();
	}
	
	double LatticeParameters[3];
	getLatticeParameters(LatticeParameters);
	
	for (int i = 0; i < N; i++)
	{
		/*	check_InBoundary.rx = ParticleNo[i].GetPosition().rx;
		if (check_InBoundary.rx != -1 * 10 * a)
		{*/
//		ParticleNo[i].positionEvolution_Constraint(evolutionType, LatticeParameters, flag_Matlab, MatlabParameters_circle );
	
		//	setNoise();
		
		ParticleNo[i].positionEvolution(evolutionType);
		ParticleNo[i].polarizationEvolution(evolutionType);
		//}
	}
}

void Lattice::printLattice()
{
	double avgPolarization=0;
	int count=1;
	cout << "Position" << "\t" << "Velocity" << "\t" << "Polarization" << endl;
	for (int i = 0; i < N; i++)
	{
		cout << ParticleNo[i].GetPosition().rx << " " << ParticleNo[i].GetPosition().ry << "\t";
		cout << ParticleNo[i].GetVelocity().vx << " " << ParticleNo[i].GetVelocity().vy << "\t";
		cout << ParticleNo[i].GetPolarization().px << " " << ParticleNo[i].GetPolarization().py <<"\t"<< ParticleNo[i].GetPolarization().modPol << endl;
	}
	//cout << "Avg Polarization - " << avgPolarization << endl;
}

void Lattice::printParameters()
{
	double Parameters[13];
	ParticleNo[0].GetParameters(Parameters); 
	ofstream out;
	out.open("Simulation Parameters.dat");
	string parameter[13];
	parameter[0] = "Timestep";
	parameter[1] = "mu1";
	parameter[2] = "mu2";
	parameter[3] = "zeta";
	parameter[4] = "k";
	parameter[5] = "v0";
	parameter[6] = "v0 at Boundary";
	parameter[7] = "Threshold v0 (if used)";
	parameter[8] = "ConfinementForce";
	parameter[9] = "Spring Stiffness when extended";
	parameter[10] = "Spring Stiffness when compressed";
	parameter[11] = "Noise strength of position";
	parameter[12] = "Noise strength of polarization";
	
	out << "Evolution Type- " << endl;
	if (evolutionType == 0){ out << "All particles are active with v0 from t=0;" << endl << endl; }
	if (evolutionType == 1){ out << "Inside particles become active only if their velocity exceeds Threshold v0" << endl << endl; }

	for (int i = 0; i < 13; i++)
	{
		cout << parameter[i] << "\t" << Parameters[i] << endl;
		out << parameter[i] << "\t" << Parameters[i] << endl;
	}
	out.close();
}

//void Particle::writeFile(){
//	double temp_force_pol;
//	temp_force_pol = Modulus(force_pol[0].fx, force_pol[0].fy);
//
//	ofstream out1;
//	out1.open("force_alongPol.dat", ios::app);
//	out1 << force_pol[0].fx << "\t" << force_pol[0].fy <<"\t"<< temp_force_pol<<endl;
//
//	ofstream out2;
//	out2.open("velocity.dat", ios::app);
//	out2 << vel.vx << "\t" << vel.vy << endl;
//
//	ofstream out3;
//	out3.open("polarization.dat", ios::app);
//	out3 << pol.px << "\t" << pol.py <<"\t"<<Modulus(pol.px,pol.py)<< endl;
//
//	ofstream out4;
//	out4.open("force.dat", ios::app);
//	out4 << force.fx << "\t" << force.fy << endl;
//
//	ofstream out5;
//	out5.open("position.dat", ios::app);
//	out5 << pos.rx << "\t" << pos.ry << endl;
//
//	ofstream out6;
//	out6.open("force_perpendicularPol.dat", ios::app);
//	out6 << force_pol[1].fx << "\t" << force_pol[1].fy<< endl;
//	
//	out1.close();
//	out2.close();
//	out3.close();
//	out4.close();
//	out5.close();
//	out6.close();
//}

void Lattice::cleanFiles(){
	ofstream out;
	out.open("lattice.dat", ios::trunc);
	out.close();
	out.open("neighbours.dat", ios::trunc);
	out.close();
	out.open("force_alongPol.dat", ios::trunc);
	out.close();
	out.open("force_perpendicularPol.dat", ios::trunc);
	out.close();
	out.open("position.dat", ios::trunc);
	out.close();
	out.open("velocity.dat", ios::trunc);
	out.close();
	out.open("polarization.dat", ios::trunc);
	out.close();
	out.open("force.dat", ios::trunc);
	out.close();
	out.open("checkEvolution.txt", ios::trunc);
	out.close();
	FILE *f1;
	f1 = fopen("checkEvolutionDouble.txt", "w");
	fclose(f1);
	FILE *f2;
	f2 = fopen("velocityDensityPlot.dat", "w");
	fclose(f2);
	/*FILE *velOrder;
	ostringstream fileName;
	fileName << "velOrder_" << EvolutionNo << ".dat";
	velOrder = fopen(fileName.str().c_str(), "a");
	velOrder = fopen("velOrder.dat", "w");
	fclose(velOrder);*/
}

void Lattice::appendFiles(int EvolutionNo){
	ofstream out;
	out.open("force_alongPol.dat", ios::app);
	out << endl;
	out.open("force_perpendicularPol.dat", ios::app);
	out << endl;
	out.close();
	out.open("position.dat", ios::app);
	out << endl;
	out.close();
	out.open("velocity.dat", ios::app);
	out << endl;
	out.close();
	out.open("polarization.dat", ios::app);
	out << endl;
	out.close();
	out.open("force.dat", ios::app);
	out << endl;
	out.close();
	FILE *f;
	f = fopen("checkEvolutionDouble.txt", "a");
	fprintf(f, "\n");
	fclose(f);
	FILE *velocityDensity;
	ostringstream fileName;
	fileName << "velocityDensityPlot" << ".dat";
	velocityDensity = fopen(fileName.str().c_str(), "a");
	fprintf(velocityDensity, "\n");
	fclose(velocityDensity);
	FILE *velOrder;
	ostringstream fileName2;
	fileName2 << "velOrder_" << EvolutionNo << ".dat";
	velOrder = fopen(fileName2.str().c_str(), "a");
	fprintf(velOrder, "\n");
	fclose(velOrder);
}

Particle Lattice::getParticle(int j)
{
	Particle P;
	P.SetPosition(ParticleNo[j].GetPosition());
	P.SetVelocity(ParticleNo[j].GetVelocity());
	P.SetPolarization(ParticleNo[j].GetPolarization().px, ParticleNo[j].GetPolarization().py);
	P.SetForce_pol(ParticleNo[j].GetForce_pol());
	return P;
}

//void Lattice::forQuiverplot(int timestep){		
//	double modForce_alongpol;
//	Force * F;
//	ostringstream fileName;
//	fileName << "quiverplot_" << timestep << ".dat";
//	ofstream out;
//	out.open(fileName.str().c_str());
//	for (int i = 0; i < N; i++)
//	{
//		F = ParticleNo[i].GetForce_pol();
//		modForce_alongpol = ParticleNo[i].Modulus(F[0].fx,F[0].fy);
//		out << ParticleNo[i].GetPosition().rx <<  "\t" << ParticleNo[i].GetPosition().ry << "\t";
//		out << ParticleNo[i].GetVelocity().vx << "\t" << ParticleNo[i].GetVelocity().vy << "\t";
//		out << ParticleNo[i].GetPolarization().px << "\t" << ParticleNo[i].GetPolarization().py <<"\t"<<ParticleNo[i].GetPolarization().modPol<< "\t";
//		out << modForce_alongpol << "\t" << F[1].fx << "\t" << F[1].fy << endl;
//	}
//	out.close();
//}

void Lattice::forQuiverplot(int timestep)
{
	double modForce_alongpol;
	Force * F;
	FILE *quiverplot;
	ostringstream fileName;
	fileName << "quiverplot_" << timestep << ".dat";
	quiverplot = fopen(fileName.str().c_str(), "w");
	Position check_InBoundary;
	for (int i = 0; i < N; i++)
	{
		check_InBoundary.rx = ParticleNo[i].GetPosition().rx;
		if (check_InBoundary.rx != -1 * 10 * a)
		{
			F = ParticleNo[i].GetForce_pol();
			modForce_alongpol = ParticleNo[i].Modulus(F[0].fx, F[0].fy);

			fprintf(quiverplot, "%.10e\t%.10e\t", ParticleNo[i].GetPosition().rx, ParticleNo[i].GetPosition().ry);
			fprintf(quiverplot, "%.10e\t%.10e\t", ParticleNo[i].GetVelocity().vx, ParticleNo[i].GetVelocity().vy);
			fprintf(quiverplot, "%.10e\t%.10e\t%.10e\t", ParticleNo[i].GetPolarization().px, ParticleNo[i].GetPolarization().py, ParticleNo[i].GetPolarization().modPol);
			fprintf(quiverplot, "%.10e\t%.10e\t%.10e\t", modForce_alongpol, F[1].fx, F[1].fy);
			fprintf(quiverplot, "%.10e\t%.10e\n", ParticleNo[i].GetForce().fx, ParticleNo[i].GetForce().fy);
		}
	}
	fclose(quiverplot);
}

void Lattice::ForceHistogramPlot(int rowNo) 
{
	int particleNo;
	Force ParticleForce;
	double ModForce;
	particleNo = L*(rowNo-1);							// L is Lattice Length
	FILE *forceHistogram;
	ostringstream fileName;
	fileName << "forceHistogramPlot_" << rowNo << ".dat";
	forceHistogram = fopen(fileName.str().c_str(), "a");
	for ( int i = particleNo; i < L*rowNo; i++) 
	{
		ParticleForce.fx = ParticleNo[i].GetForce().fx;
		ParticleForce.fy = ParticleNo[i].GetForce().fy;
		ModForce = ParticleNo[i].Modulus(ParticleForce.fx, ParticleForce.fy);
		fprintf(forceHistogram, "%.10e\t%.10e\t%.10e\n", ParticleForce.fx, ParticleForce.fy, ModForce);
	}
	fclose(forceHistogram);
}

void Lattice::newForceHistogram(int rowNo, int timestep)
{
	int particleNo;
	Force ParticleForce;
	double ModForce;
	particleNo = L*(rowNo - 1);							// L is Lattice Length
	FILE *forceHistogram;
	ostringstream fileName;
	fileName << "forceHistogramPlot_" << rowNo <<"_"<<timestep<<".dat";
	forceHistogram = fopen(fileName.str().c_str(), "a");
	for (int i = particleNo; i < L*rowNo; i++)
	{
		ParticleForce.fx = ParticleNo[i].GetForce().fx;
		ParticleForce.fy = ParticleNo[i].GetForce().fy;
		ModForce = ParticleNo[i].Modulus(ParticleForce.fx, ParticleForce.fy);
		fprintf(forceHistogram, "%.10e\t%.10e\t%.10e\n", ParticleForce.fx, ParticleForce.fy, ModForce);
	}
	fclose(forceHistogram);
}

void Lattice::VelOrderParameter(int EvolutionNo)
{
	Force * F;
	FILE *velOrder;
	Velocity currentVel;
	Position currentPos;
	int i;
	int n;
	double x0;	//centre of circle
	double y0;  //centre of circle
	double diff_x;
	double diff_y;
	double cos_theta;
	double sin_theta;
	double diff_modulus;
	double VelOrder=0;
	double temp_VelOrder;
	double VelModulus;
	x0 = MatlabParameters_circle[0];
	y0 = MatlabParameters_circle[1];
	for (int n = 0; n < N; n++)
	{
		currentPos = ParticleNo[n].GetPosition();
		currentVel = ParticleNo[n].GetVelocity();
		diff_x = (currentPos.rx - x0);
		diff_y = (currentPos.ry - y0);
		if (diff_x == 0 || diff_y == 0){ diff_modulus = 1; }
		else{ diff_modulus = ParticleNo[n].Modulus(diff_x, diff_y); }
		cos_theta = diff_x / diff_modulus;
		sin_theta = diff_y / diff_modulus;
		//Rotating the cos_theta and sin_theta vector by 90 in clockwise direction gives [sin_theta -cos_theta]
		temp_VelOrder = ParticleNo[n].DotProduct(currentVel.vx, currentVel.vy, sin_theta, -1.0 * cos_theta);
		VelModulus = ParticleNo[n].Modulus(currentVel.vx,currentVel.vy);
		if (VelModulus==0){	VelModulus = 1.0;}
		temp_VelOrder = temp_VelOrder / VelModulus;
		VelOrder = VelOrder + temp_VelOrder;
		//cout << n << " " << VelOrder << endl;
	}
	ostringstream fileName;
	fileName << "velOrder_" << EvolutionNo << ".dat";
	velOrder = fopen(fileName.str().c_str(), "a");
	fprintf(velOrder, "%.10e\t", VelOrder);
	fclose(velOrder);
}

void Lattice::VelocityContourPlot(int iterationNo)			// plot the avg vel of all the particles in a row
{
	double totalModVelocity;
	double avgModVelocity;
	Velocity avgVelocity;
	Velocity totalVelocity;
	int rowNo = 1;
	int i = 0;
	int z;
	while (i<N)
		{
			if (i%L==0)				// L is total no. of particles in a row
			{
				totalModVelocity = 0;
				totalVelocity.vx = 0;
				totalVelocity.vy = 0;
				z = i + L;			// to count the particles that are only in a single row
				while (i<z)
				{
					totalVelocity.vx = totalVelocity.vx + ParticleNo[i].GetVelocity().vx;
					totalVelocity.vy = totalVelocity.vy + ParticleNo[i].GetVelocity().vy;
					totalModVelocity = totalModVelocity + ParticleNo[i].GetVelocity().modVel;
					i++;
				}
				avgVelocity.vx = totalVelocity.vx / L;
				avgVelocity.vy = totalVelocity.vy / L;
				avgModVelocity = totalModVelocity / L;
				
				FILE *velocityDensity;
				ostringstream fileName;
				fileName << "velocityDensityPlot" << ".dat";
				velocityDensity = fopen(fileName.str().c_str(), "a");
				fprintf(velocityDensity, "%d\t%d\t%.10e%\t%.10e\t%.10e\n", iterationNo, rowNo, avgVelocity.vx, avgVelocity.vy, avgModVelocity);
				fclose(velocityDensity);
				
				rowNo++;
			}
		}
}		

double Lattice::ran2(long *idum)
{
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		idum2 = (*idum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ1;
			*idum = IA1*(*idum - k*IQ1) - k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ1;
	*idum = IA1*(*idum - k*IQ1) - k*IR1;
	if (*idum < 0) *idum += IM1;
	k = idum2 / IQ2;
	idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j = iy / NDIV;
	iy = iv[j] - idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp = AM*iy) > RNMX) return RNMX;
	else return temp;
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
}

double Lattice::gasdev(long *idum)	//Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum) as the source of uniform deviates.
{
	//double ran2(long *idum);
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;
	if (*idum < 0) iset = 0; //Reinitialize.
	if (iset == 0) {	//We don’t have an extra deviate handy, so
		do {
			v1 = 2.0*ran2(idum) - 1.0; //pick two uniform numbers in the square exv2 =
			v2 = 2.0*ran2(idum) - 1.0; //tending from - 1 to + 1 in each direction,
			rsq = v1*v1 + v2*v2; //see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); //and if they are not, try again.
		fac = sqrt(-2.0*log(rsq) / rsq);
		//Now make the Box - Muller transformation to get two normal deviates.Return one and save the other for next time.
		gset = v1*fac;
		iset = 1; //Set flag.
		return v2*fac;
	}
	else {		//We have an extra deviate handy,
		iset = 0; //so unset the flag,
		return gset; //and return it.
	}
}

void Lattice::test_checkNeighbours(int i, Position * neighbours) {
	ofstream out;
	out.open("neighbours.dat", ios::app);
	out << ParticleNo[i].GetPosition().rx << "  " << ParticleNo[i].GetPosition().ry << " ";
	for (int j = 0; j < 6; j++)
	{
		out << neighbours[j].rx - ParticleNo[i].GetPosition().rx << "  " << neighbours[j].ry - ParticleNo[i].GetPosition().ry << " ";
	}
	out << endl;
	out.close();
}

void Lattice::test_checkNeighbours_ConstraintHex()
{
	int TotalNeighbours;
	if (flag_Matlab==0)
	{
		TotalNeighbours = 6;
	}
	else if (flag_Matlab==1)
	{
		TotalNeighbours = MaxNeighbours_Matlab;
	}
	ofstream out;
	out.open("neighbours_Constraint.dat");
	Position current;
	for (int i = 0; i < N; i++)
	{
		current = ParticleNo[i].GetPosition();
		out << ParticleNo[i].GetPosition().rx << "  " << ParticleNo[i].GetPosition().ry << " ";
		for (int j = 0; j < TotalNeighbours; j++)
		{
			out << ParticleNo[i].NeighbourArray[j].pos.rx - ParticleNo[i].GetPosition().rx << "  " << ParticleNo[i].NeighbourArray[j].pos.ry - ParticleNo[i].GetPosition().ry << " ";
			//if (j>3){ cout << i << "  " <<j<<"  "<< ParticleNo[i].NeighbourArray[j].pos.rx << endl; }
		}
		out << endl;
	}
	out.close();
}

void Lattice::test_calculateForce(int i, Position * neighbours)
{
	if (i == 13)
	{
		printf("to check Force \n");
		printf("%.10e\t%.10e\n", ParticleNo[i].GetPosition().rx, ParticleNo[i].GetPosition().ry);
		for (int j = 0; j < 6; j++)
		{
			printf("%.10e\t%.10e\n", neighbours[j].rx, neighbours[j].ry);
		}
		printf("%.10e\t%.10e\n", ParticleNo[i].GetForce().fx, ParticleNo[i].GetForce().fy);
	}

}

void Particle::print_checkEvolution(int particleNo, int time_iteration, int count)
{
	double temp_force_pol;
	temp_force_pol = Modulus(force_pol[0].fx, force_pol[0].fy);  // modulus of force acting along polarizatoin
	FILE *f;
	f = fopen("checkEvolutionDouble.txt", "a");
	//	ofstream out;
	//	out.open("checkEvolution.txt", ios::app);
	if (count == 0)
	{
		//		out << "Time Step = " << h << ",  " << "Spring Constant = " << k << ",  " << "V0 = " << v0 << ",  " << "Zeta = " << zeta << ",  " << "Mu1 = " << mu1 << ",  " << "Mu2 = " << mu2 << endl;
		//		out << endl;
		fprintf(f, "Time Step = %.1e \t Spring Constant = %2.2f \t V0 = %2.2f \t Zeta = %2.2f \t Mu1 = %2.2f \t Mu2 = %2.2f \n \n", h, k, v0, zeta, mu1, mu2);
	}
	//	out << "Particle No. " << particleNo << ",  " << "At time iteration " << time_iteration << endl;
	fprintf(f, "Particle No. %d At time iteration %d \n", particleNo, time_iteration);
	//	out << "Position"<<"\t"<< pos.rx << "\t" << pos.ry << endl;
	fprintf(f, "Position\t%.16e\t%.10e\n", pos.rx, pos.ry);
	//	out << "Velocity"<<"\t"<<vel.vx << "\t" << vel.vy << endl;
	//	out << "Polarization and modP" << "\t" << pol.px << "\t" << pol.py << "\t" << Modulus(pol.px, pol.py) << endl;
	fprintf(f, "Polarization and modP\t%.16e\t%.10e\t%.10e\n", pol.px, pol.py, Modulus(pol.px, pol.py));
	//	out << "Force along Polarization and its magnitude"<<"\t"<<force_pol[0].fx << "\t" << force_pol[0].fy <<"\t" <<temp_force_pol << endl;
	fprintf(f, "Magnitude of Force along Polarization\t%.16e\n", temp_force_pol);
	//	out << "Force perpendicular to Polarization" << "\t"<<force_pol[1].fx << "\t" << force_pol[1].fy << endl;
	fprintf(f, "Force perpendicular to Polarization\t%.16e\t%.10e\n", force_pol[1].fx, force_pol[1].fy);
	//	out.close();
	fclose(f);
}

//void Lattice::Evolve_ConstraintHex()
//{
//	getNeighbours_ConstraintHex();
//	for (int i = 0; i < N; i++)
//	{
//		calculateForce_ConstraintHex(i);
//		ParticleNo[i].calculateForce_pol();
//	}
//	for (int i = 0; i < N; i++)
//	{
//		/*	check_InBoundary.rx = ParticleNo[i].GetPosition().rx;
//		if (check_InBoundary.rx != -1 * 10 * a)
//		{*/
//		ParticleNo[i].positionEvolution(evolutionType);
//		ParticleNo[i].polarizationEvolution(evolutionType);
//		//}
//	}
//}

//Position Particle::positionEvolution_Constraint(int evolutionType, double * LatticeParameters, int flag_Matlab, double * MatlabParameters_Circle)
//{
//	double force_alongPol;		// modulus of component of force parallel to polarization
//	double modVel;				// magnitude of velocity
//	double temp_v0;
//	double sin60 = sqrt(3) / 2;
//	double cos60;
//	cos60 = 0.5;
//	double L, B;
////	double position;           // current position of particle from centre
//	double radius_inner;			   // radius of confinement/ECM	
//	double radius_outer;			   // radius of confinement/ECM	
//	double x0;
//	double y0;
//	double a0;
//	double diff_x;
//	double diff_y;
//	double cos_theta;
//	double sin_theta;
//	double diff_modulus;
//	double extension;
//	double flag_BoundaryCross= 0;	 // to check whether particle has crossed the boundary
//	//	double temp_cout;
//	L = LatticeParameters[0];
//	B = LatticeParameters[1];
//	//cout << "Lattice Parameter: L is " << L << " B is " << B << endl;
//	force_alongPol = DotProduct(pol.px, pol.py, force.fx, force.fy);
//	modVel = Modulus(vel.vx, vel.vy);
//	if (flag_Matlab==0)
//	{
//		x0 = L - 1;
//		y0 = (L - 1)*sin60;
//		radius_inner = L - 1;
//		radius_outer = L + 1;
//	}
//	else if (flag_Matlab==1)
//	{
//		x0 = MatlabParameters_Circle[0];
//		y0 = MatlabParameters_Circle[1];
//		radius_inner = MatlabParameters_Circle[2];
//		radius_outer = MatlabParameters_Circle[3];
//	}
//	a0 = radius_outer - radius_inner;
//	diff_x = (pos.rx - x0);
//	diff_y = (pos.ry - y0);
//	if (diff_x == 0 || diff_y == 0){	diff_modulus = 1;	}
//	else{ diff_modulus = Modulus(diff_x, diff_y); }
//	cos_theta = diff_x / diff_modulus;
//	sin_theta = diff_y / diff_modulus;
//	//position = diff_modulus*diff_modulus;
//	extension = diff_modulus - radius_inner;
//	if (extension>0)
//	{
//		flag_BoundaryCross = 1.0;
//		/*cout << "In confinement force condition" << endl;
//		cin>>temp_cout;*/
//	}
//	if (evolutionType == 0)
//	{
//		vel.vx = (mu1*force_alongPol + v0)*pol.px + mu2*force_pol[1].fx - flag_BoundaryCross*ConfinementForce*extension*cos_theta;			//diff_x/diff_modulus gives the direction to which force is added
//		vel.vy = (mu1*force_alongPol + v0)*pol.py + mu2*force_pol[1].fy - flag_BoundaryCross*ConfinementForce*extension*sin_theta;
//		pos.rx = pos.rx + h*vel.vx;
//		pos.ry = pos.ry + h*vel.vy;
//	}
//	if (evolutionType == 1)
//	{
//		if (modVel >= 0.25)
//		{
//			count = 1;
//		}
//		if (count == 0){ temp_v0 = 0; }
//		else{ temp_v0 = v0; }
//		vel.vx = (mu1*force_alongPol + temp_v0)*pol.px + mu2*force_pol[1].fx - flag_BoundaryCross*ConfinementForce*extension*cos_theta;
//		vel.vy = (mu1*force_alongPol + temp_v0)*pol.py + mu2*force_pol[1].fy - flag_BoundaryCross*ConfinementForce*extension*sin_theta;
//		pos.rx = pos.rx + h*vel.vx;
//		pos.ry = pos.ry + h*vel.vy;
//	}
//	return pos;
//}

//void Lattice::forQuiverplot1(){
//	ofstream out;
//	out.open("quiverplot1.dat");
//	for (int i = 0; i < N; i++)
//	{
//		out << ParticleNo[i].GetPosition().rx << "  " << ParticleNo[i].GetPosition().ry << " ";
//		out << ParticleNo[i].GetVelocity().vx << "  " << ParticleNo[i].GetVelocity().vy << endl;
//	}
//}
//void Lattice::forQuiverplot2(){
//	ofstream out;
//	out.open("quiverplot2.dat");
//	for (int i = 0; i < N; i++)
//	{
//		out << ParticleNo[i].GetPosition().rx << "  " << ParticleNo[i].GetPosition().ry << " ";
//		out << ParticleNo[i].GetVelocity().vx << "  " << ParticleNo[i].GetVelocity().vy << endl;
//	}
//}

//Velocity * Particle::calculateVelocity_pol()
//{
//	double temp_vx_pol = 0, temp_vy_pol = 0;		//temp_vx is parallel to Polarization and temp_vy is perpendicular
//	double temp_vx = 0, temp_vy = 0;
//	Velocity temp;
//
//	temp_vx_pol = DotProduct(pol.px, pol.py, vel.vx, vel.vy);	//mod of component of velocity along polarization (mod(Pol)=1)
//	velocity_pol[0].vx = temp_vx_pol*pol.px;		//component of velocity along to polarization- in x-direction
//	velocity_pol[0].vy = temp_vx_pol*pol.py;		//component of velocity along polarization- in y-direction
//	temp_vy_pol = Modulus((temp_vx - vel.vx), (temp_vy - vel.vy));  //mod of component of velocity perpendicular to polarization (mod(Pol)=1) 
//	velocity_pol[1].vx = vel.vx - temp_vx_pol*pol.px;
//	velocity_pol[1].vy = vel.vy - temp_vx_pol*pol.py;
//	temp.vx = velocity_pol[0].vx + velocity_pol[1].vx;
//	temp.vy = velocity_pol[0].vy + velocity_pol[1].vy;
//	SetVelocity(temp);
//	
//	ofstream out1;
//	out1.open("vel_test.txt");
//	out1 << temp_vx_pol<<endl;
//	out1 << temp.vx << " " << temp.vy << endl;
//	out1 << pol.px << " " << pol.py << endl;
//	return GetVelocity_pol();
//}
