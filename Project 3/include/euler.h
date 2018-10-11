#ifndef EULER_H
#define EULER_H

class Euler
{
    private:
        double dt;

    public:
        Euler (double dt);
        void integrateOneStep(class SolarSystem &system);
};

#endif
