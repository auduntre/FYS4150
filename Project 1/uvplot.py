import matplotlib.pyplot as plt
import numpy as np

def uvplot():
    ns = [10, 100, 1000]

    for n in ns:
        u = np.loadtxt("u" + str(n) + ".txt")
        v = np.loadtxt("v" + str(n) + ".txt")

        x = np.linspace(0, 1, n+2)

        # TODO: add legends!
        plt.plot(x[1:-1], u)
        plt.plot(x[1:-1], v, 'r--')
        plt.show()


if __name__ == "__main__":
    uvplot()
