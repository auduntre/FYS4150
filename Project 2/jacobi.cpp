#include <armadillo>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>


arma::uword maxoffdiag (arma::mat X)
{
    arma::mat Y = arma::eye<arma::mat>(X.n_rows, X.n_cols);

    // Set diag to zero and find index to extreme value in matrix
    Y = arma::abs(X - (X % Y)); // % == Schur product 
    return Y.index_max();
}



TEST_CASE ("Implementation of maximum off-diagonal absolute value method", 
           "[maxoffdiag]") 
{
    arma::mat A;
    int n = 100;
    double big = -123.123;
    
    A.randu(n, n);
    A(n/3, n/4) = big; // set maximum off-diagonal value
    
    SECTION ("Testing if right index") {    
        REQUIRE (maxoffdiag(A) == n*(n/4) + n/3); // Check if right index
    }
    SECTION ("Testing if right value") {
        REQUIRE (A(maxoffdiag(A)) == big);
    }
}
