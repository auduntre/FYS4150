#include <iostream>
#include <cstdlib>
#include <cmath>
#include <armadillo>
#include <string>

int main(int argc, char **argv)
{
    arma::vec poisson_solver (arma::vec a, arma::vec b, arma::vec c, 
                              arma::vec v, unsigned int n);
    double u_exact (double x);

    arma::vec n_vec(3); 
    n_vec << 10 << 100 << 1000;


    std::string vfile = "v_";
    std::string ufile = "u_";
    std::string refile = "re_";
    std::string ext = ".txt";
    std::string tmp;

    std::string strarr[3] = {vfile, ufile, refile}; 
    int n;
 
    for (int i = 0; i < n_vec.n_elem; i++) {
        n = n_vec[i];
        tmp = std::to_string(n);

        arma::vec a(n);
        arma::vec b(n);
        arma::vec c(n);
        arma::vec v(n);

        a.fill(-1.0);
        b.fill(2.0);
        c.fill(-1.0);

        v = poisson_solver(a, b, c, v, n);

        arma::vec x = arma::linspace<arma::vec> (0.0, 1.0, n+2);
        arma::vec u(n);
        arma::vec re(n);

        for (int j = 0; j < n; j++) {
            u[j] = u_exact(x[j+1]);
            re[j] = std::fabs((v[j] - u[j]) / u[j]);
        }

        for (std::string str: strarr) {
            str.append(tmp);
            str.append(ext);
        }

        v.save(vfile, arma::raw_ascii);
        u.save(ufile, arma::raw_ascii);
        re.save(refile, arma::raw_ascii); 
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
