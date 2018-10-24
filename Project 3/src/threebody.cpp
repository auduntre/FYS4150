#include <cstdlib>

#include "solarsystem.h"
#include "euler.h"
#include "verlet.h"

int main (int argc, char **argv)
{
    void threebody (double massj);

    double massjs[4] = {0.95e-3, 0.95e-2, 0.95e-1, 0.95};

    for (double mj: massjs) {
        threebody(mj);
    }

    return 0;
}


void threebody (double massj)
{
    SolarSystem sol;

    arma::vec3 pos = {0, 0, 0};
    arma::vec3 vel = {0, -0.007231252772102223, 0};
    double mass = 1.0;

    CelestialBody &sun = sol.createCelestialBody(pos, vel, mass);

    pos = {5.20, 0.0, 0.0};
    vel = {0.0, 7.592003385453352, 0.0};
    mass  = massj;

    CelestialBody &jupiter = sol.createCelestialBody(pos, vel, mass);

    pos = {1, 0.0, 0.0};
    vel = {0.0, 2 * M_PI, 0.0};
    mass  = 3e-6;

    CelestialBody &earth = sol.createCelestialBody(pos, vel, mass);

    int nTimesteps = 1000000;
    double endTime = 10.0;
    double dt = endTime / nTimesteps;

    Verlet integrator(dt);
    integrator.initialStep(sol);

    arma::mat posMat(3, nTimesteps, arma::fill::zeros); 
    posMat.col(0) = earth.position;

    for (int i = 1; i < nTimesteps; i++) {
        integrator.integrateOneStep(sol);

        posMat.col(i) = earth.position - sun.position;
    }

    posMat.save("../positions/positions" + std::to_string(massj) + ".txt",
                arma::raw_ascii);

    earth.position.print();

    return;
}
