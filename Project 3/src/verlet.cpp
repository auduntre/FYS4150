#include "verlet.h"
#include "solarsystem.h"

Verlet::Verlet (double deltaT): dt(deltaT)
{
}


void Verlet::integrateOneStep (SolarSystem &system)
{
    system.calculateForcesAndEnergy();

    for (CelestialBody &body : system.bodies()) {
        arma::vec3 fwrdPos = 2.0 * body.position - body.past_position
                           + (body.force / body.mass) * (this->dt * this->dt);


        body.velocity = (fwrdPos - body.past_position) / (2 * this->dt);
        body.past_position = body.position;
        body.position = fwrdPos;
    }
}
