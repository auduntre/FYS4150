#include <iostream>
#include <cmath>
#include <cstdlib>

#include "solarsystem.h"
#include "euler.h"


int main (int argc, char **argv)
{
    int nTimesteps = 1000;

    if (argc >= 2) {
        nTimesteps = std::atoi(argv[1]);
    }

    SolarSystem sol;

    arma::vec3 pos = {0, 0, 0};
    arma::vec3 vel = {0, 0, 0};
    double mass = 1.0;

    CelestialBody &sun = sol.createCelestialBody(pos, vel, mass);

    pos = {1, 0, 0};
    vel = {0, 2*M_PI, 0};
    mass  = 3e-6;

    CelestialBody &earth = sol.createCelestialBody(pos, vel, mass);

    std::vector<CelestialBody> &bodies = sol.bodies();
}