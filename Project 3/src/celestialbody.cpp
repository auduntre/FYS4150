#include "celestialbody.h"

CelestialBody::CelestialBody (arma::vec3 pos, arma::vec3 vel, double m) :
    position(pos),
    velocity(vel),
    mass(m)
{
}


CelestialBody::CelestialBody (double x, double y, double z, double vx,
                              double vy, double vz, double m) : 
    mass(m)
{
    this->position = {x, y, z};
    this->velocity = {vx, vy, vz};
    this->past_position = this->position;
}


void CelestialBody::resetForce ()
{
    this->force = {0, 0, 0};
}
