#include "solarsystem.h"
#include <iostream>

SolarSystem::SolarSystem () : 
    ke(0),
    pe(0)
{
}


CelestialBody &SolarSystem::createCelestialBody (arma::vec3 position, 
                                                 arma::vec3 velocity,
                                                 double mass)
{
    this->bods.push_back(CelestialBody(position, velocity, mass));
    return this->bods.back();
}


void SolarSystem::calculateForcesAndEnergy (bool forceAndPE, bool angMomAndKE)
{
    if (angMomAndKE) {
        this->ke = 0;
        this->angMom = {0, 0, 0};
    }
    
    // Reset forces on all bodies
    if (forceAndPE) {
        this->pe = 0;
        
        for (CelestialBody &body : this->bods) {
            body.past_force = body.force;
            body.resetForce();
        }
    }
    
    for (int i = 0; i < this->numberOfBodies(); i++) {
        CelestialBody &body1 = this->bods[i];

        if (forceAndPE) {
            for (int j = i+1; j < this->numberOfBodies(); j++) {
                CelestialBody &body2 = this->bods[j];
                arma::vec3 deltaRVector = body1.position - body2.position;
                arma::vec3 normRVector = arma::normalise(deltaRVector, 2);

                // Euclidean distance is the L2-norm
                double dr = arma::norm(deltaRVector, 2);

                // Calculate potential energy
                double U = - GravConst * body1.mass * body2.mass / dr;
                this->pe += U;

                // Calculate the force
                double drPow1 = this->forcePower(dr, body1);
                double drPow2 = this->forcePower(dr, body2);
                body1.force += (U / drPow1) * normRVector;
                body2.force -= (U / drPow2) * normRVector;
            }
        }
    
        if (angMomAndKE) {
            // Calculate kinetic enerfy
            this->ke += 0.5 * body1.mass 
                      * arma::dot(body1.velocity, body1.velocity);
        
            // Calculate angular momentum
            body1.computeOrbitalAngMom();
            this->angMom += body1.mass * body1.orbAngMom;

        }
    }
}


double SolarSystem::forcePower (double r, CelestialBody body)
{
    return r;
}


int SolarSystem::numberOfBodies () const
{
    return this->bods.size();
}


double SolarSystem::totalEnergy () const
{
    return this->ke + this->pe;
}


double SolarSystem::potentialEnergy () const
{
    return this->pe;
}


double SolarSystem::kineticEnergy () const
{
    return this->ke;
}


arma::vec3 SolarSystem::angularMomentum () const
{
    return this->angMom;
}


std::vector<CelestialBody> &SolarSystem::bodies ()
{
    return this->bods;
}


void SolarSystem::writeToFile (std::string filename)
{
    if (this->ofile.good()) {
        this->ofile.open(filename.c_str(), std::ofstream::out);
        
        if(!this->ofile.good()) {
            std::cout << "Error opening file " << filename << ". Aborting!" << std::endl;
            std::terminate();
        }
    }

    this->ofile << this->numberOfBodies() << std::endl;
    this->ofile << "Position of bodies below:" << std::endl;

    for(CelestialBody &body : this->bods) {
        this->ofile << "1 " << body.position(0) << " " 
                    << body.position(1) << " " 
                    << body.position(2) << "\n";
    }

    this->ofile.close();
}
