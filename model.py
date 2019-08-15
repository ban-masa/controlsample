import numpy as np
import matplotlib.pyplot as plt

m1 = 1.0
m2 = 1.0
k = 100.0
dt = 0.01
d = 1.0

A = np.array([[1.0, dt, 0.0, 0.0],
              [-k * dt / m1, 1.0, k * dt / m1, 0.0],
              [0.0, 0.0, 1.0, dt],
              [k * dt / m2, 0.0, -k * dt / m2, 1.0 - d * dt / m2]])

B = np.array([[0.0, dt / m1, 0.0, 0.0]]).T

c = np.array([[0, 0, 1, 0]])

def calc_next(x, u):
    global A
    global B
    return np.dot(A, x) + B * u

def test1():
    timelist = []
    datalist = []
    x = np.array([[0, 0, 0, 0]]).T
    for i in range(10000):
        timelist.append(i * dt)
        u = -1 * (x[2][0] - 5.0)
        x = calc_next(x, u)
        datalist.append(x[2][0])
    plt.plot(timelist, datalist, label="test1")
    plt.legend()
    plt.show()

def riccatiIteration(A, B, Q, R):

def calc_gain():


def test2():
    timelist = []
    datalist = []
    x = np.array([[0,0,0,0]]).T
    K = calc_gain()
    for i in range(10000):
        timelist.append(i * dt)
        u = np.dot(K, x)[0][0]
        x = calc_next(x, u)
        datalist.append(x[2][0])
    plt.plot(timelist, datalist, label="test2")
    plt.legend()
    plt.show()

    

test1()
