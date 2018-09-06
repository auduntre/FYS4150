#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def uv_plot():
    ns = [10, 100, 1000]

    for n in ns:
        u = np.loadtxt("results/u" + str(n) + ".txt")
        v = np.loadtxt("results/v" + str(n) + ".txt")

        x = np.linspace(0, 1, n+2)

        # TODO: add legends!
        plt.plot(x[1:-1], u)
        plt.plot(x[1:-1], v, 'r--')
        plt.show()

def re_table():
    ns = [10 ** i for i in range(1, 8)]
    hs = [1.0 / (n + 1) for n in ns]
    res = np.loadtxt("results/re_max.txt")

    print("{:9}  |  {:14}".format("N:", "Re_max:"))
    print("-----------|----------------")

    row_form = "{:>9}  |  {:>14}"
    for n, re in zip(ns, res):
        print(row_form.format(n, re))


if __name__ == "__main__":
    uv_plot()
    re_table()
