#ifndef HAROSCILL_H
#define HAROSCILL_H

#include <armadillo>
#include <cmath>

typedef unsigned uint;

class Harmonic_oscillator
{
    private:
        uint N;

        double rhoN;
        double h;
        double dconst;
        double econst;
        double eps;

        arma::vec rho;
        arma::vec Vi;

    public:
        Harmonic_oscillator(double rhoN, uint N);

        arma::vec one_electron ();
        void two_electrons ();

        arma::vec analy_eigval ();
};

#endif