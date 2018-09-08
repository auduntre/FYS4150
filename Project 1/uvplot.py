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
        plt.savefig("results/uvplot" + str(n) + ".png")
        plt.figure()


def re_table():
    res = np.loadtxt("results/re_max.txt")
    hs = np.loadtxt("results/h.txt");
    ns = [10**i for i in range(1, len(hs) + 1)]

    with open("results/re_table.txt", "w") as re_table:
        header = "{:9}  |  {:17}  |  {:14}\n"
        re_table.write(header.format("N:", "h:", "log(Re_max):"))
        re_table.write("-----------|--------------------|----------------\n")

        row_form = "{:>9}  |  {:.17g}  |  {:>14}\n"
        for n, re, hs in zip(ns, res, hs):
            re_table.write(row_form.format(n, hs, re))


if __name__ == "__main__":
    uv_plot()
    re_table()
