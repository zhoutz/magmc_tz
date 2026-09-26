import matplotlib.pyplot as plt
import numpy as np

settings = [
    ("output/ft07_beta075_twist10_E.txt", "-"),
    ("output/ft07_beta075_twist10_O.txt", "--"),
]


def f(fname, ls):
    data = np.loadtxt(fname)
    muk = data[:, 1]
    plt.hist(muk, bins=100,  histtype="step", label=fname, ls=ls)


for fname, ls in settings:
    f(fname, ls)

plt.show()
