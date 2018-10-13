#include "solarsystem.h"
#include <iostream>

SolarSystem::SolarSystem () : 
    ke(0),
    pe(0),
    beta(2.0)
{
}


CelestialBody &SolarSystem::createCelestialBody (arma::vec3 position, 
                                                 arma::vec3 velocity,
                                                 double mass)
{
    this->bods.push_back(CelestialBody(position, velocity, mass));
    return this->bods.back();
}


void SolarSystem::calculateForcesAndEnergy ()
{
    this->calculateForcesAndPE();
    this->calculateAngMomAndKE();
}


void SolarSystem::calculateForcesAndPE ()
{
    this->pe = 0;
    
    // Reset forces on all bodies
    for (CelestialBody &body : this->bods) {
        body.past_force = body.force;
        body.resetForce();
    }
    
    for (int i = 0; i < this->numberOfBodies(); i++) {
        CelestialBody &body1 = this->bods[i];

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
            body1.force += (U / dr) * normRVector;
            body2.force -= (U / dr) * normRVector;
        }
    }
}


void SolarSystem::calculateAngMomAndKE ()
{
    this->ke = 0;
    this->angMom = {0, 0, 0};
    
    for (int i = 0; i < this->numberOfBodies(); i++) {
        CelestialBody &body1 = this->bods[i];

        // Calculate kinetic enerfy
        this->ke += 0.5 * body1.mass 
                  * arma::dot(body1.velocity, body1.velocity);
        
        // Calculate angular momentum
        this->angMom += arma::cross(body1.position,
                                    body1.mass * body1.velocity);
    }
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
}


arma::vec3 SolarSystem::angularMomentum () const
{
    return this->angMom;
}


std::vector<CelestialBody> &SolarSystem::bodies ()
{
    return this->bods;
}
