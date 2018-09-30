#include <cmath>
#include <iostream>
#include <iomanip>

#include "jacobi.h"
#include "toeplitz.h"


constexpr double pi() {
    return 4 * std::atan(1);
}


int main(int argc, char **argv)
{
    arma::vec analy_eigval(double a, double d, uint N);

    double eps = 1.0E-16;
    double a;
    double d;

    int N = 100;


}


bool obtain_analy()


arma::vec analy_eigval(double a, double d, uint N)
{
    arma::vec js = arma::linspace<arma::vec>(1, N-1, N-1);
    arma::vec lambda(N-1);

    lambda = d + 2 * a * arma::cos(js * pi() / (N + 1));

    return lambda;
}