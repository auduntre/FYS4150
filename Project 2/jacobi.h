#ifndef JACOBI_H
#define JACOBI_H

#include <armadillo>

typedef unsigned uint;

arma::mat jacobi (arma::mat A, double eps, uint maxiter);

arma::mat rotate (arma::mat B, uint k, uint l);

double maxoffdiag (arma::mat X);

#endif
