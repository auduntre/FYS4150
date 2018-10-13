#ifndef EULER_H
#define EULER_H


class Euler
{
    private:
        double dt; // == delta_t == dt

    public:
        Euler (double deltaT);
        
        void integrateOneStep (class SolarSystem &system);
        void integrateNtimes (class SolarSystem &system, int nTimes);

        // get dt
        double h() const;
};

#endif