#include "verlet.h"
#include "solarsystem.h"

Verlet::Verlet (double deltaT): Euler(deltaT)
{
}


void Verlet::integrateOneStep (SolarSystem &system)
{
    system.calculateForcesAndEnergy();

    for (CelestialBody &body : system.bodies()) {
        arma::vec3 fwrdPos = 2.0 * body.position - body.past_position
                           + (body.force / body.mass) * this->h() * this->h();


        body.velocity = (fwrdPos - body.past_position) / (2 * this->h());
        body.past_position = body.position;
        body.position = fwrdPos;
    }
}


void Verlet::initialStep(class SolarSystem &system)
{
    for (CelestialBody &body : system.bodies()) {
        body.past_position = body.position;
    }

    Euler::integrateOneStep(system);
}


void Verlet::integrateNtimes(class SolarSystem &system, int nTimes)
{
    this->initialStep(system);

    for (int i = 1; i < nTimes; i++) {
        this->integrateOneStep(system);
    }
}