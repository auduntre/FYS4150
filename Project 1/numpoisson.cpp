#include <iostream>
#include <cstdlib>
#include <cmath>
#include <armadillo>
#include <string>


int main(int argc, char **argv)
{
    arma::vec poisson_solver (arma::vec f_vec, arma::vec a, arma::vec b,
                              arma::vec c, unsigned int n);
    double u_exact (double x);
    double f (double x); 

    arma::vec n_vec(3);
    n_vec << 10 << 100 << 1000;
    int n;

    for (int i = 0; i < n_vec.n_elem; i++) {
        n = n_vec[i];

        arma::vec x = arma::linspace<arma::vec> (0.0, 1.0, n+2);
        
        arma::vec f_vec(n);
        arma::vec a(n);
        arma::vec b(n);
        arma::vec c(n);
        
        a.fill(-1.0);
        b.fill(2.0);
        c.fill(-1.0);
        
        arma::vec v(n);
        arma::vec u(n);
        arma::vec re(n);

        for (int j = 0; j < n; j++) {
            u[j] = u_exact(x[j+1]);
            f_vec[j] = f(x[j+1]); 
        }

        v = poisson_solver(f_vec, a, b, c, n);

        for (int j = 0; j < n; j++) {
            re[j] = std::fabs((v[j] - u[j]) / u[j]);
        }

        if (n == 100) {
            for (int k = 0; k < n; k++) {
                std::cout << "re = " << re[k] <<
                          ", u[" << k << "] = " << u[k] <<
                          ", v[" << k << "] = " << v[k] <<
                          std::endl;
            }
            std::cout << "DONE" << std::endl;
        }
    }
    return 0;
}



arma::vec poisson_solver(arma::vec f_vec, arma::vec a, arma::vec b,
                         arma::vec c, unsigned int n)
{
    arma::vec x = arma::linspace<arma::vec>(0, 1.0, n+2);
    arma::vec v(n);

    double h = 1.0 / ((double) n + 1.0); // or (x[1] - x[0]);
    double hh = h * h;

    //Forward subst
    f_vec[0] = f_vec[0] * hh;
    for (int i = 1; i < n; i++) {
        b[i] = b[i] - (a[i-1] * c[i-1]) / b[i-1];
        f_vec[i] = (f_vec[i] * hh) - (a[i-1] * f_vec[i-1]) / b[i-1];
    }

    //Backward subst
    v[n-1] = f_vec[n-1] / b[n-1];
    for (int i = n-2; i >= 0; i--) {
        v[i] = (f_vec[i] - c[i] * v[i+1]) / b[i];
    }

    return v;
}


double f(double x)
{
    return 100.0 * std::exp(-10.0 * x);

}


double u_exact(double x)
{
    return 1.0 - (1.0 - std::exp(-10.0)) * x - std::exp(-10 * x);
}
