#include "solarsystem.h"
#include <iostream>

SolarSystem::SolarSytem() : 
    ke(0),
    pe(0)
{
}


CelestialBody &SolarSystem::createCelestialBody(vec3 position, vec3 velocity,
                                                double mass)
{
    this->bods.push_back(CelestialBody(poistion, velocity, mass));
    return this->bods.back();
}


void SolarSystem::calculateForcesAndEnergy ()
{
    this->ke = 0;
    this->pe = 0;
    this->angMom = {0, 0, 0};

    // Reset forces on all bodies
    for (CelestialBody &body : this->bods) {
        body.force = {0, 0, 0};
    }

    for (int i = 0; i < this->numberOfBodies(); i++) {
        CelestialBody &body1 = this->bods[i];

        for (int j = i+1; j < this->numberOfBodies; j++) {
            CelestialBody &body2 = this->bods[j];
            vec3 deltaRVector = body1.position - body2.position;
            
            // Euclidean distance is the L2-norm
            double dr = arma::norm(deltaRVector, 2);
            
            // Calculate the force and potential energy here
        }

        this->ke += 0.5 * body1.mass 
                  * arma::dot(body1.velocity, body1.velocity);
}
