#include <iostream>
#include <math.h>
#include "Particle.h"
using namespace std;

#define PI 3.14159265358979323846

Particle::Particle()
{
	h = 0.0001;		//delta-t (time interval) between 2 updates in Euler method
	mu1 = 15.0;		
	mu2 = 1.0;		
	zeta = 5.0;
	v0 = 5.0;
	v0_Boundary = 7.5;
	k = 3.0;
	k_extension=5.5;
	k_compression = 15.0;   //spring constant when spring is compressed
	thresholdVel = 0.1;
	D[0] = 0;	//Noise strength for position
	D[1] = 0;  //Noise strength for polarization
	NoisePreFactor[0] = D[0] * sqrt(h);  //Noise pre-factor for position
	NoisePreFactor[1] = D[1] * sqrt(h);  //Noise pre-factor for polarization
	noisePol = 0;
	count[2] = 0;	// to check if the threshold has been crossed for motility
	boundary = 0;   // will become 1 if particle is at Boundary
	NeighbourArray = new Neighbours[8];
	NeighbourNo_Matlab = new int[8];
	ConfinementForce = k;
	//Neighbours();
}

//Neighbours::Neighbours()
//{
//	EquilibriumDistance = 1;
//}

Position::Position()
{
	a = 1.0;
	rx = -1*10 * a;
	ry = -1*10 * a;
}

Velocity::Velocity()
{
	vx = 0.0;
	vy = 0.0;
	modVel = 0.0;
}

Polarization::Polarization()
{
	px = 0.0;
	py = 0.0;
	modPol = 1.0;
}


Force::Force()
{
	fx = 0.0;
	fy = 0.0;
}

Noise::Noise()
{
	nx = 0.0;
	ny = 0.0;
}

void Particle::GetParameters(double * Parameters)
{
	Parameters[0] = h;
	Parameters[1] = mu1;
	Parameters[2] = mu2;
	Parameters[3] = zeta;
	Parameters[4] = k;
	Parameters[5] = v0;
	Parameters[6] = v0_Boundary;
	Parameters[7] = thresholdVel;
	Parameters[8] = ConfinementForce;
	Parameters[9] = k_extension;
	Parameters[10] = k_compression;
	Parameters[11] = D[0];
	Parameters[12] = D[1];
}

//double Particle::Getv0_Boundary(){
//	return v0_Boundary;
//}

void Particle::SetPosition(Position p){
	pos.rx = p.rx;
	pos.ry = p.ry;
}

void Particle::SetPosition(double x, double y){
	pos.rx = x;
	pos.ry = y;
}

Position Particle::GetPosition(){
	return pos;
}

Velocity Particle::GetVelocity(){
	return vel;
}

Polarization Particle::GetPolarization(){
	return pol;
}

void Particle::SetVelocity(Velocity v)
{
	vel.vx = v.vx;
	vel.vy = v.vy;
}

void Particle::SetPolarization(double x, double y){
	pol.px = x;
	pol.py = y;
}

void Particle::SetForce(double x, double y){
	force.fx = x;
	force.fy = y;
}

Force Particle::GetForce(){
	return force;
}

//void Particle::SetForce_pol(double x, double y){
//	force_pol.fx = x;
//	force_pol.fy = y;
//}

Force * Particle::GetForce_pol(){
	return force_pol;
}

void Particle::SetForce_pol(Force * Force_pol) {
	force_pol[0].fx = Force_pol[0].fx;
	force_pol[0].fy = Force_pol[0].fy;
	force_pol[1].fx = Force_pol[1].fx;
	force_pol[1].fy = Force_pol[1].fy;
}

void Particle::SetNoise(Noise * temp_noise){
	noise[0].nx = temp_noise[0].nx;
	noise[0].ny = temp_noise[0].ny;
	noise[1].nx = temp_noise[1].nx;
	noise[1].ny = temp_noise[1].ny;
}

Velocity * Particle::GetVelocity_pol(){
	return velocity_pol;
}


void Particle::setNeighbourNo_Matlab(int NeighbourNo, int j)
{
	NeighbourNo_Matlab[j] = NeighbourNo;
}
int Particle::getNeighbourNo_Matlab(int j)
{
	return NeighbourNo_Matlab[j];
}
//void Particle::setNeighbours()
//{
//
//}

double Particle::DotProduct(double x1, double y1, double x2, double y2){
	double dotproduct;
	dotproduct = x1*x2 + y1*y2;
	return dotproduct;
}

double Particle::Modulus(double x, double y){
	double mod;
	mod = sqrt(x*x + y*y);
	return mod;
}

//
//int main(){
//	Particle P;
//	Position p;
//	p.x = 1;
//	p.y = 2;
//	P.SetPosition(p);
//	P.SetPosition(2, 3);
//	P.GetPosition();
//	cin >> p.x;
//	return 0;
//}