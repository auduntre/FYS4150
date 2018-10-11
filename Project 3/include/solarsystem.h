#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "celestialbody.h"

#include <fstream>
#include <vector>
#include <string>


class SolarSystem
{
    private:
        std::vector<CelestialBody> m_bodies;
        std::ofstream m_file;
        
        vec3 m_angularMomentum;
        double m_kineticEnergy;
        double m_potentialEnergy;

    public:
        SolarSystem();

        CelestialBody &createCelestialBody (vec3 position, vec3 velocity,
                                            double mass);
        void calculateForcesAndEnergy ();
        int numberOfBodies () const;

        double totalEnergy () const;
        double potentialEnergy () const;
        double kineticEnergy () const;
        void writeToFile (std::string filename);
        arma::vec angularMomentum () const;
        std::vector<CelestialBody> &bodies ();
};

#endif
