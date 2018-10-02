#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np


def plot_two_electron():
    filenamevec = "results/eigvector100rhoN5omega"
    filenameval = "results/eigvalue100rhoN5omega"

    ext = ".txt"
    omegas = ["0.010000", "0.500000", "1.000000", "5.000000"]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    first = True
    for omega in omegas:
        ev = np.loadtxt(filenamevec + omega + ext)
        el = np.loadtxt(filenameval + omega + ext)
        psi2 = el * ev**2

        if first:
            x = np.linspace(0, 5.0, psi2.size)
            psi2inter = psi2
            first = False

        cl = r"$\omega_r$ = " + omega
        plt.plot(x, psi2, label=cl)
    
    plt.title(r"Relative energy of ground state")
    plt.ylabel(r"$|\psi (\rho)|^2$")
    plt.xlabel(r"$\rho$")
    
    plt.legend()
    plt.show()


    plt.


if __name__ == '__main__':
    plot_two_electron()