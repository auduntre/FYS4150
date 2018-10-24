#include <cstdlib>
#include <string>

#include "betasolarsystem.h"
#include "euler.h"
#include "verlet.h"


int main(int argc, char **argv)
{
    void betaRun (double beta);

    arma::vec betas = {2.00, 2.50, 2.75, 2.90, 2.95, 3.00};
    
    for (double beta: betas) {
        betaRun(beta);
    }

    return 0;
}


void betaRun(double beta)
{

    BetaSolarSystem sol(beta);

    double mass_sun = 1.0;
    arma::vec3 pos = {0.0, 0.0, 0.0};
    arma::vec3 vel = {0.0, -1.8849555921538758e-05, 0.0};

    CelestialBody &sun = sol.createCelestialBody(pos, vel, mass_sun);
    
    double mass_earth = 3.0e-6;
    pos = {1.0, 0.0, 0.0};
    vel = {0.0, 6.28, 0.0};
    
    CelestialBody &earth = sol.createCelestialBody(pos, vel, mass_earth);
    
    int nTimesteps = 500000;
    double endTime = 5.0;
    double dt = endTime / nTimesteps;

    Verlet integrator(dt);

    arma::mat posMat(3, nTimesteps, arma::fill::zeros);
    posMat.col(0) = earth.position;
    
    for (int i = 1; i < nTimesteps; i++) {
        integrator.integrateOneStep(sol);
        posMat.col(i) = earth.position;
    }

    posMat.save("../positions/beta" + std::to_string(beta) + ".txt",
                arma::raw_ascii);
    
    return;
}
