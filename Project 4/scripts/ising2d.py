# Code for the two-dimensional Ising model with periodic boundary conditions
# taken from lecture notes in FYS3150 

from matplotlib import pyplot as plt
from numba import jit
import numpy as np
import sys, math

@jit(nopython=True)
def periodic (i, limit, add):
    """
    Choose correct matrix index with periodic
    boundary conditions

    Input:
    - i:     Base index
    - limit: Highest \"legal\" index
    - add:   Number to add or subtract from i
    """
    return (i + limit + add) % limit


@jit(nopython=True)
def monteCarlo(temp, size, trials):
    """
    Calculate the enerInvalid usage of BoundFunction(array.item for array(int64, 2d, C)) with parameters (int64, int64)
    * parameterizegy and magnetization
    (\"straight\" and squared) for a given temperature

    Input:
    - temp:   Temperature to calculate for
    - size:   dimension of square matrix
    - trials: Monte-carlo trials (how many times do we
                                  flip the matrix?)

    Output:
    - E_av:       Energy of matrix averaged over trials, normalized to spins**2
    - E_variance: Variance of energy, same normalization * temp**2
    """

    #Setup spin matrix, initialize to ground state
    spin_matrix = np.zeros( (size,size), np.int8) + 1

    #Create and initialize variables
    E = 0
    E_av = E2_av = 0

    #Setup array for possible energy changes
    w = np.zeros(17, np.float64)
    for de in range(-8, 9, 4): #include +8
        w[de+8] = math.exp(-de/temp)

    #Calculate initial energy
    for j in range(size):
        for i in range(size):
            E -= spin_matrix[i, j] * \
                 (spin_matrix[periodic(i, size, -1), j] + spin_matrix[i, periodic(j, size, 1)])

    #Start metropolis MonteCarlo computation
    for i in range(trials):
        #Metropolis
        #Loop over all spins, pick a random spin each time
        for s in range(size**2):
            x = int(np.random.random()*size)
            y = int(np.random.random()*size)
            deltaE = 2 * spin_matrix[x,y] * \
                     (spin_matrix[periodic(x, size, -1), y] +\
                      spin_matrix[periodic(x, size, 1),  y] +\
                      spin_matrix[x, periodic(y, size, -1)] +\
                      spin_matrix[x, periodic(y, size, 1)])
            if np.random.random() <= w[deltaE+8]:
                #Accept!
                spin_matrix[x,y] *= -1
                E += deltaE

        #Update expectation values
        E_av    += E
        E2_av   += E**2

    E_av       /= float(trials);
    E2_av      /= float(trials);
    #Calculate variance and normalize to per-point and temp
    E_variance  = (E2_av-E_av*E_av)/float(size*size*temp*temp);
    #Normalize returned averages to per-point
    E_av       /= float(size*size);

    return (E_av, E_variance)


def main():
    """Main program"""

    # values of the lattice, number of Monte Carlo cycles and temperature domain
    size        = 20
    trials      = 100000
    temp_init   = 1.8
    temp_end    = 2.6
    temp_step   = 0.1


    temps = np.arange(temp_init, temp_end+temp_step/2, temp_step, float)
    Dim = np.size(temps)
    energy = np.zeros(Dim)
    heatcapacity = np.zeros(Dim)
    temperature = np.zeros(Dim)
    
    for i, temp in enumerate(temps):
        (E_av, E_variance) = monteCarlo(temp,size,trials)
        temperature[i] = temp
        energy[i] = E_av
        heatcapacity[i] = E_variance
    
    plt.figure(1)
    plt.subplot(211)
    plt.axis([1.8,2.6,-2.0, -1.0])
    plt.xlabel(r'Temperature $J/(k_B)$')
    plt.ylabel(r'Average energy per spin  $E/N$')
    plt.plot(temperature, energy, 'b-')
    plt.subplot(212)
    plt.axis([1.8,2.6, 0.0, 2.0])
    plt.plot(temperature, heatcapacity, 'r-')
    plt.xlabel(r'Temperature $J/(k_B)$')
    plt.ylabel(r'Heat capacity per spin  $C_V/N$')
    plt.savefig('energycv.pdf')
    #plt.show()


if __name__ == "__main__":
    main()
