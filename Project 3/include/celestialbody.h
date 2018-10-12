#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H

#include <armadillo>

// Fixed vectors in 3 dimensional space

class CelestialBody
{
    public:
        arma::vec3 position;
        arma::vec3 velocity;
        arma::vec3 force;
        double mass;

        CelestialBody (arma::vec3 pos, arma::vec3 vel, double m);
        CelestialBody (double x, double y, double z, double vx, double vy,
                       double vz, double mass);
        
        void resetForce();
};

#endif
