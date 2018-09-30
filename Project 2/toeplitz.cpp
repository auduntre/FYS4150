#include "toeplitz.h"

arma::mat sym_tridiag(double a, double d, uint N)
{
    arma::vec diags = arma::zeros<arma::vec>(N);
    diags(0) = d;
    diags(1) = a;

    return arma::toeplitz(diags);
}