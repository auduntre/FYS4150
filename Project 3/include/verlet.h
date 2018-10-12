#ifndef VERLET_H
#define VERLET_H

class Verlet
{
    private:
        double dt;

    public:
        Verlet (double deltaT);

        void integrateOneStep(class SolarSystem &system);
};

#endif