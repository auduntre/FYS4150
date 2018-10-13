#include "euler.h"
#include "solarsystem.h"


Euler::Euler (double deltaT) : dt(deltaT) 
{
}


void Euler::integrateOneStep (SolarSystem &system)
{
    for (CelestialBody &body : system.bodies()) {
        body.position += body.velocity * this->dt;
        body.velocity += (body.force / body.mass) * this->dt;
    }

    system.calculateForcesAndEnergy();
}


void Euler::initialStep(class SolarSystem &system)
{
    system.calculateForcesAndEnergy();
    
    this->integrateOneStep (system);

}

void Euler::integrateNtimes (class SolarSystem &system, int nTimes)
{
    this->initialStep(system);

    for (int i = 1; i < nTimes; i++) {
        this->integrateOneStep(system);
    }
}


double Euler::h() const
{
    return this->dt;
}
