#ifndef HAROSCILL_H
#define HAROSCILL_H

#include <armadillo>
#include <cmath>

typedef unsigned uint;


/*
 * Class for studying Harmonic Oscilliation
 */
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
        /*
         * Setting up the harmoinc oscillator.
         *
         * @param rhoN: maximum value of the dimensionless variable
         * rho = r / alhpa for the relative wave function psi(rho).
         * @param N: Number of resolution points for the study.
         */
        Harmonic_oscillator(double rhoN, uint N);

        /*
         * Study of the case with one electron.
         */
        arma::vec one_electron ();
        
        /*
         * Study of the case with two electrons.
         *
         * @param *Ev: eigenvector Armadillo-matrix with dimensions 
         * (N-1)x(N-1). Must be submitted as an indentity matrix.
         * @param omega: Strength of the oscillator potential.
         * @param inter: If there it is a interacting system or not.
         * @return Eigenvalues of the system in an Armadillo vector.
         */  
        arma::vec two_electrons (arma::mat *Ev, double omega, bool inter);

        /*
         * Returns the analytical eigenvalues for the one electron case.
         *
         * @return Eigenvalues in an Armadillo vector.
         */
        arma::vec analy_eigval_one_election();
};

#endif
