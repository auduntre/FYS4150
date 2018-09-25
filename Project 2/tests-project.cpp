#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "jacobi.h"

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
