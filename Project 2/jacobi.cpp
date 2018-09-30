#include <iostream>
#include <cmath>

#include "jacobi.h"


arma::mat jacobi (arma::mat A, double eps, uint maxiter) 
{
    void maxoffdiag (arma::mat X, uint *i, uint *j);

    arma::mat B = A;
    arma::mat S(A.n_rows, A.n_rows);

    uint maxi;
    uint maxj;
    uint iter = 0;

    double tau;
    double root_tau;
    double tan_theta;
    double cos_theta;
    double sin_theta;

    maxoffdiag(B, &maxi, &maxj);

    while (std::fabs(B(maxi, maxj)) > eps && maxiter > iter) {
        tau = (B(maxj, maxj) - B(maxi, maxi)) / (2.0 * B(maxi, maxj));

        root_tau = std::sqrt(1 + (tau * tau));

        if (tau > 0) {
            tan_theta = 1.0 / (tau + root_tau);
        } else {
            tan_theta = -1.0 / (-tau + root_tau);
        }
        
        cos_theta = 1.0 / root_tau;
        sin_theta = tan_theta * cos_theta;

        S.eye();
        S(maxi, maxj) = - sin_theta;
        S(maxj, maxi) = sin_theta;
        S(maxi, maxi) = cos_theta;
        S(maxj, maxj) = cos_theta;

        B = S.t() * B * S;

        maxoffdiag(B, &maxi, &maxj);
        iter++;
    }

    return B;
}


void maxoffdiag (arma::mat X, uint *i, uint *j)
{
    // Set diag to zero and find index of extreme value in matrix
    uint ij = arma::abs(arma::trimatl(X) - arma::diagmat(X)).index_max();

    *j = ij / X.n_rows;
    *i = ij - (*j) * X.n_rows;

    return;
}