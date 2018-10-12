#ifndef VERLET_H
#define VERLET_H

#include "euler.h"

class Verlet : public Euler
{

    public:
        Verlet (double deltaT);

        void integrateOneStep(class SolarSystem &system);
        void initialStep (class SolarSystem &system);
        void integrateNtimes (class SolarSystem &system, int nTimes);
};

#endif