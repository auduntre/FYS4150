#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H

#include <armadillo>

// Fixed vectors in 3 dimensional space
typedef arma::vec::fixed<3> vec3;

class CelestialBody
{
    public:
        vec3 position;
        vec3 velocity(3);
        vec3 force(3);
        double mass;

        CelestialBody (vec3 position, vec3 velocity, double mass);
        CelestialBody (double x, double y, double z, double vx, double vy,
                       double vz, double mass);
        
        void resetForce();
};

#endif
