#include "catch.hpp"
#include <iostream>

#include "jacobi.h"


bool equalcolumns(arma::mat A, arma::mat B, double eps)
{
    bool column_k;

    for (uint k = 0; k < A.n_cols; k++) {
        column_k = false;

        for (uint l = 0; l < B.n_cols; l++) {
            if (arma::approx_equal(A.col(k), B.col(l), "absdiff", eps)) {
                column_k = true;
            }
            else if (arma::approx_equal(-A.col(k), B.col(l), "absdiff", eps))
            {
                column_k = true;
            }
            
        }

        if (!column_k) {
            return false;
        }
    }
    return true;
}


bool orthogonal(arma::mat R, double eps)
{
    arma::mat Reye(R.n_rows, R.n_cols);
    Reye.eye(); 

    if (arma::approx_equal(R.t() * R, Reye, "absdiff", eps)) {
        return true;
    } else {
        return false;
    }
}


TEST_CASE ("Implementation of maximum off-diagonal absolute value method", 
           "[maxoffdiag]") 
{
    arma::mat A;
    uint n = 1000;
    uint i;
    uint j;
    double big = (double) n * (double) n;
    
    // Make A symmetric
    A.randu(n, n);
    A = A.t() * A;

    A(n/4, n/3) = big; // set maximum off-diagonal value some place;
    
    SECTION ("Testing if right index") {
        maxoffdiag(A, &i, &j);
        REQUIRE ((i == n/4 && j == n/3)); // Check if right index
    }
    SECTION ("Testing if right value") {
        REQUIRE (maxoffdiag(A, &i, &j) == big);
    }
}


TEST_CASE ("Test implementation of Jacobi's rotation algorithm", 
           "[jacobi]")
{
    bool equalcolumns(arma::mat A, arma::mat B, double eps);
    bool orthogonal(arma::mat R, double eps);

    double eps = 1.0E-12;
    double dps = 1.0E-8;

    int n = 100;

    arma::mat A(n ,n);
    arma::mat Ev(n, n);

    Ev.eye();
    A.randu();
    A = A.t() * A;

    arma::mat eigvec;
    arma::vec eigval;
    arma::vec jacobival;
    
    arma::eig_sym(eigval, eigvec, A);
    jacobival = arma::sort(jacobi(A, &Ev, eps));

    SECTION ("Testing if right eigenvalues") {
        REQUIRE (arma::approx_equal(jacobival, eigval, "absdiff", dps));
    }
    SECTION ("Testing if right eigenvetors") {
        REQUIRE (equalcolumns(Ev, eigvec, dps));
    }
    SECTION ("Testing if eigenvectors have orthogonality preserved") {
        REQUIRE (orthogonal(Ev, dps));
    }
}