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
    ns = np.logspace(1, np.size(res), np.size(res), dtype=np.int)

    with open("results/re_table.txt", "w") as re_table:
        re_table.write("{:9}  |  {:14}\n".format("N:", "log(Re_max):"))
        re_table.write("-----------|----------------\n")

        row_form = "{:>9}  |  {:>14}\n"
        for n, re in zip(ns, res):
            re_table.write(row_form.format(n, re))


if __name__ == "__main__":
    uv_plot()
    re_table()
