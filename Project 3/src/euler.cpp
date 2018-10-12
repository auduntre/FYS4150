#include "euler.h"
#include "solarsystem.h"

Euler::Euler (double deltaT) : dt(deltaT) 
{
}


void Euler::integrateOneStep (SolarSystem &system)
{
    system.calculateForcesAndEnergy();

    for (CelestialBody &body : system.bodies()) {
        body.position += body.velocity * this->dt;
        body.velocity += (body.force / body.mass) * this->dt;
    }
}
