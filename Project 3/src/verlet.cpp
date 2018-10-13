#include "verlet.h"
#include "solarsystem.h"


Verlet::Verlet (double deltaT): Euler(deltaT)
{
    this->dtdt = this->h() * this->h(); 
}


void Verlet::integrateOneStep (SolarSystem &system)
{
    system.calculateForcesAndEnergy();
    
    for (CelestialBody &body : system.bodies()) {
        body.position = body.position + body.velocity * this->h()
                      + 0.5 * (body.force / body.mass) * this->dtdt;
    }
    
    system.calculateForcesAndEnergy();

    for (CelestialBody &body : system.bodies()) {
        body.velocity = body.velocity + 0.5 * this->h()
                      * (body.past_force / body.mass 
                      + body.force / body.mass);
    }
}


void Verlet::initialStep(class SolarSystem &system)
{
    system.calculateForcesAndEnergy();
    
    for (CelestialBody &body : system.bodies()) {
        body.position = body.position + body.velocity * this->h()
                      + 0.5 * (body.force / body.mass) * this->dtdt;
    }
    
    system.calculateForcesAndEnergy();

    for (CelestialBody &body : system.bodies()) {
        body.velocity = body.velocity + 0.5 * this->h()
                      * (body.past_force / body.mass 
                      + body.force / body.mass);
    }
}


void Verlet::integrateNtimes(class SolarSystem &system, int nTimes)
{
    this->initialStep(system);

    for (int i = 1; i < nTimes; i++) {
        this->integrateOneStep(system);
    }
}

