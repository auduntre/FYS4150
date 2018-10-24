#include <cstdlib>
#include <string>

#include "solarsystem.h"
#include "euler.h"
#include "verlet.h"


int main(int argc, char **argv)
{
    void betaRun (double beta);

    for (double beta: betas) {  
    }

    return 0;
}


void betaRun(double beta)
{
    arma::vec betas = {2.00, 2.50, 2.75, 2.90, 2.95, 3.00};

    double mass_earth = 3.0e-6;
    arma::vec3 pos_earth = {1.0, 0.0, 0.0};
    arma::vec3 vel_earth = {0.0, 6.28, 0.0};
    
    CelestialBody &sun = sol.createCelestialBody(pos, vel, mass);
    
    double mass_sun = 1.0;
    arma::vec3 pos_sun = {0.0, 0.0, 0.0};
    arma::vec3 vel_sun = - vel_earth * mass_earth / mass_sun;

    CelestialBody &earth = sol.createCelestialBody(pos, vel, mass);

    double nTimesteps = 5000000;
    double endTime = 5.0;

    return;
}
