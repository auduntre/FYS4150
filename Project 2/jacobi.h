#ifndef JACOBI_H
#define JACOBI_H

#include <armadillo>
#include <cmath>

typedef unsigned uint;

arma::vec jacobi (arma::mat A, arma::mat *Ev, double eps);

arma::mat rotate (arma::mat B, arma::mat *Ev, uint k, uint l);

double maxoffdiag (arma::mat X, uint *i, uint *j);

#endif