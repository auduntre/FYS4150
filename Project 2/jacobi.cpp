#include <armadillo>
#include <cmath>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>


arma::uword maxoffdiag (arma::mat X)
{
    // Set diag to zero and find index of extreme value in matrix
    return arma::abs(arma::trimatu(X)).index_max();
}



TEST_CASE ("Implementation of maximum off-diagonal absolute value method", 
           "[maxoffdiag]") 
{
    arma::mat A;
    int n = 1000;
    double big = n * n;
    
    // Make A symmetric
    A.randu(n, n);
    A = A.t() * A;

    A(n/4, n/3) = big; // set maximum off-diagonal value some place;
    
    SECTION ("Testing if right index") {    
        REQUIRE (maxoffdiag(A) == n*(n/3) + n/4); // Check if right index
    }
    SECTION ("Testing if right value") {
        REQUIRE (A(maxoffdiag(A)) == big);
    }
}
