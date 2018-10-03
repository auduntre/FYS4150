#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np


def plot_two_electron():
    """Plots the results from the two electron case."""
    filenamevec = "results/eigvector100rhoN5omega"
    filenameval = "results/eigvalue100rhoN5omega"

    ext = ".txt"
    omegas = ["0.010000", "0.500000", "1.000000", "5.000000"]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    ##### COMAPRING DIFFRENET OSCILLATOR POTENTIALS #####
    first = True
    for omega in omegas:
        ev = np.loadtxt(filenamevec + omega + ext)
        psi2 = ev**2

        if first:
            rho = np.linspace(0, 5.0, psi2.size)
            psi2inter = psi2
            first = False

        cl = r"$\omega_r$ = " + omega
        plt.plot(rho, psi2, label=cl)
    
    plt.title(r"Relative energy of ground state")
    plt.ylabel(r"$|\psi (\rho)|^2$")
    plt.xlabel(r"$\rho$")
    
    plt.legend()
    plt.show()

    ##### COMPARING INTERACTING AND NON-INTERACTING CAS #####
    ev = np.loadtxt("results/noninteractionvector.txt")
    psi2non = ev**2

    plt.title(r"With and without Coulomb interaction $(\omega = 0.01)$")
    plt.ylabel(r"$|\psi (\rho)|^2$")
    plt.xlabel(r"$\rho$")

    plt.plot(rho, psi2inter, label="interacting")
    plt.plot(rho, psi2non, label="non-interacting")

    plt.legend()
    plt.show()


if __name__ == '__main__':
    plot_two_electron()
