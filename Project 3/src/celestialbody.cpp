#include "celestialbody.h"

CelestialBody::CelestialBody (arma::vec3 pos, arma::vec3 vel, double m) :
    position(pos),
    velocity(vel),
    mass(m)
{
    this->resetForce(); // Make sure forces is initialised
}


CelestialBody::CelestialBody (double x, double y, double z, double vx,
                              double vy, double vz, double m) : 
    mass(m)
{
    this->position = {x, y, z};
    this->velocity = {vx, vy, vz};
}


void CelestialBody::resetForce ()
{
    this->force = {0, 0, 0};
}
