# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:34:34 2024

@author: HFKJ059
"""

import numpy as np
from matplotlib import pyplot as plt


def read_1d(filename):
    with open("./data/" + filename + ".txt", "r") as file:
        data = file.read()
    data = data.split()
    for i in range(len(data)):
        data[i] = float(data[i])
    data = np.array(data)
    return data

def read_2d(filename):
    data = []
    with open("./data/" + filename + ".txt", "r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            data_temp = line.split()
            data.append(data_temp)
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = float(data[i][j])
    data = np.array(data)
    return data

def read(index):
    x = read_1d("x")
    y = read_1d("y")
    phi = read_2d("phi" + index)
    return x, y, phi

def plot(x, y, phi, title):
    X, Y = np.meshgrid(x, y)
    plt.contourf(X, Y, phi.T, levels=1000, cmap="jet")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.colorbar()
    plt.show()

def run(index, title):
    x, y, phi = read(index)
    plot(x, y, phi, title)

def main():
    run("1", "explicit Euler")
    run("2", "implicit Euler")
    run("3", "Crank-Nicolson")
    run("4", "three-time-level")


main()
