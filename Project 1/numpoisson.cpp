#include <iostream>
#include <cstdlib>
#include <cmath>
#include <armadillo>


int main(int argc, char **argv)
{
    arma::vec poisson_solver (arma::vec a, arma::vec b, arma::vec c, 
                              arma::vec v, unsigned int n);
    double u_exact (double x);

    arma::vec n_vec(3); 
    n_vec << 10 << 100 << 1000;

    double n;
    
    for (int i = 0; i < n_vec.n_elem; i++) {
        n = n_vec[i];

        arma::vec a(n);
        arma::vec b(n);
        arma::vec c(n);
        arma::vec v(n);

        a.fill(-1.0);
        b.fill(2.0);
        c.fill(-1.0);

        v = poisson_solver(a, b, c, v, n);
        
        if (n == 100) {
            arma::vec x = arma::linspace<arma::vec> (0, 1, n+2);
            arma::vec u(n);
            for (int k = 1; k < n; k++) {
                u[k] = u_exact(x[k]);

                double re = std::fabs((v[k-1] - u[k]) / u[k]);

                std::cout << "u(" << x[k] << ") = " << u[k] <<
                    ", v[" << k-1 << "] = " << v[k-1] <<
                    ", re = " << re << std::endl;
            }
        }
    }

    return 0;
}


arma::vec poisson_solver(arma::vec a, arma::vec b, arma::vec c, 
                         arma::vec v, unsigned int n)
{
    double f (double x);

    arma::vec x = arma::linspace<arma::vec>(0, 1.0, n+2);
    arma::vec fh_vec(n);

    double h = x[1] - x[0];

    //Forward subst
    for (int i = 1; i < n; i++) {
        b[i] = b[i] - (a[i-1] * c[i-1]) / b[i-1];
        fh_vec[i] = (f(x[i+1]) * h * h) - (a[i-1] * fh_vec[i-1]) / b[i-1];
    }

    //Backward subst
    v[n-1] = fh_vec[n-1] / b[n-1];
    for (int i = n-2; i >= 0; i--) {
        v[i] = (fh_vec[i] - c[i]*v[i+1]) / b[i];
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
