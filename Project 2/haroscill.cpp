#include "haroscill.h"
#include "jacobi.h"
#include "toeplitz.h"

#include <iostream>

Harmonic_oscillator::Harmonic_oscillator(double rhoN, uint N)
{
    this->N = N;
    this->rhoN = rhoN;
    
    this->h = rhoN / N;
    this->dconst = 2.0 / (this->h * this->h);
    this->econst = -1.0 / (this->h * this->h);

    this->eps = 1E-8;

    this->rho = arma::linspace<arma::vec>(h, rhoN - h, N-1);
    this->Vi = (this->rho % this->rho);
}


arma::vec Harmonic_oscillator::one_electron()
{
    uint N = this->N;
    
    arma::vec eigval = this->analy_eigval_one_election();
    arma::vec jeigval(N-1);
    arma::vec tmp(4);

    arma::mat A(N-1, N-1);
    arma::mat Ev = arma::eye<arma::mat>(N-1, N-1);

    Tridiag td(this->econst, this->dconst, N-1);
    td.set_variable_diags(this->Vi);
    
    A = td.get_matrix();
    jeigval = arma::sort(jacobi(A, &Ev, eps));

    for (uint i = 0; i < tmp.n_elem; i++) {
        tmp(i) = std::fabs(jeigval(i) - eigval(i));
    }

    return tmp;
}


arma::vec Harmonic_oscillator::two_electrons(arma::mat *Ev, double omega,
                                            bool inter)
{
    uint N = this->N;

    arma::mat A(N-1, N-1);

    arma::vec jeigval(N-1);
    arma::vec dvar = omega * this->Vi;

    if (inter == true) {
        dvar = dvar + 1 / this->rho;
    }

    Tridiag td(this->econst, this->dconst, N-1);
    td.set_variable_diags(dvar);

    A = td.get_matrix();
    jeigval = jacobi(A, Ev, eps);

    return jeigval;
}


arma::vec Harmonic_oscillator::analy_eigval_one_election()
{
    return arma::linspace<arma::vec>(3.0, 3.0 + 4.0*(this->N-2), this->N-1);
}