import numpy as np
import matplotlib.pyplot as plt

m1 = 2.0
m2 = 2.0
k = 10.0
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
    return timelist, datalist

def solveRiccatiArimotoPotter(A, B, Q, R):
    dim_x = A.shape[0]
    dim_u = B.shape[1]
    Ham = np.zeros((2 * dim_x, 2 * dim_x))
    Ham[0:dim_x,0:dim_x] = A
    Ham[0:dim_x, dim_x:2*dim_x] = np.dot(np.dot(-B, np.linalg.inv(R)), B.T)
    Ham[dim_x:2*dim_x, 0:dim_x] = -Q
    Ham[dim_x:2*dim_x, dim_x:2*dim_x] = -A.T
    w, v = np.linalg.eig(Ham)

    eigvec = np.zeros((2 * dim_x, dim_x), dtype=np.complex)
    j = 0
    for i in range(2 * dim_x):
        if w[i].real < 0.0:
            eigvec[:,j] = v[:,i]
            j = j+1

    Vs_1 = eigvec[0:dim_x, 0:dim_x]
    Vs_2 = eigvec[dim_x:2*dim_x, 0:dim_x]
    P = np.dot(Vs_2, np.linalg.inv(Vs_1)).real
    return P

def solveRiccatiIteration(A, B, Q, R):
    dim_x = A.shape[0]
    P = np.zeros((dim_x, dim_x))
    for i in range(10000):
        R_btPb_inv = np.linalg.inv(R + np.dot(np.dot(B.T, P), B))
        tmp_pa = np.dot(P, A)
        K = np.dot(np.dot(R_btPb_inv, B.T), tmp_pa)
        prev_P = np.dot(A.T, tmp_pa) + Q - np.dot(np.dot(np.dot(A.T, P), B), K)
        if (np.abs(P - prev_P) < 5.0e-10).all():
            print "ok"
            return P
        P = prev_P
    return P


def calc_gain(A, B, Q, R, Time):
    P = solveRiccatiIteration(A, B, Q, R)
    R_btPb_inv = np.linalg.inv(R + np.dot(np.dot(B.T, P), B))
    tmp_PA = np.dot(P, A)
    K = np.dot(np.dot(R_btPb_inv, B.T), tmp_PA)
    A_minusBKT = (A - np.dot(B, K)).T
    f_list = []
    dim_x = A.shape[0]
    temp =np.dot(R_btPb_inv, B.T)
    for i in range(int(Time / dt)):
        fi = np.dot(temp, Q)
        temp = np.dot(temp, A_minusBKT)
        f_list.append(fi)

    return K, f_list

def calc_u(K, f_list, x, ref_x):
    u = np.dot(-K, x)
    for f in f_list:
        u = u + np.dot(f, ref_x)
    return u

def test2():
    dim_x = A.shape[0]
    dim_u = B.shape[1]
    #Q = np.zeros((dim_x, dim_x))
    Q = np.identity(dim_x)
    Q[0,0] = 0.1
    Q[1,1] = 0.1
    Q[2,2] = 100.0
    Q[3,3] = 0.1
    R = np.zeros((dim_u, dim_u))
    R[0,0] = 1
    K, f_list = calc_gain(A, B, Q, R, 2.0)

    timelist = []
    datalist = []
    ref_x = np.array([[0, 0, 5, 0]]).T
    x = np.array([[0,0,0,0]]).T
    for i in range(10000):
        timelist.append(i * dt)
        u = calc_u(K, f_list, x, ref_x)
        x = calc_next(x, u[0,0])
        datalist.append(x[2][0])
    return timelist, datalist
    
def view_result():
    timelist1, datalist1 = test1()
    timelist2, datalist2 = test2()
    plt.plot(timelist1, datalist1, label="test1")
    plt.plot(timelist2, datalist2, label="test2")
    plt.legend()
    plt.show()

view_result()
