#include "celestialbody.h"

CelestialBody::CelestialBody (arma::vec3 pos, arma::vec3 vel, double m) :
    position(pos),
    velocity(vel),
    mass(m)
{
    this->resetForce(); // Make sure forces is initialised
    this->computeOrbitalAngMom();
}


CelestialBody::CelestialBody (double x, double y, double z, double vx,
                              double vy, double vz, double m) : 
    mass(m)
{
    this->position = {x, y, z};
    this->velocity = {vx, vy, vz};
    
    this->resetForce(); // Make sure forces is initialised
    this->computeOrbitalAngMom();
}


void CelestialBody::computeOrbitalAngMom ()
{
    // Per unit mass
    this->orbAngMom = arma::cross(this->position, this->velocity);
}


void CelestialBody::resetForce ()
{
    this->force = {0, 0, 0};
}
