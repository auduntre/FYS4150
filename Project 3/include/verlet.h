#ifndef VERLET_H
#define VERLET_H

#include "euler.h"


class Verlet : public Euler
{
    private:
        double dtdt;

    public:
        Verlet (double deltaT);

        void integrateOneStep(class SolarSystem &system);
};

#endif
