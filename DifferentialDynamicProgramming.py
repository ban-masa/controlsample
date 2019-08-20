import numpy as np

#Qx = lx + fx^T * Vx(i+1)
#Qu = lu + fu^T * Vx(i+1)
#Qxx = lxx + fx^T * Vxx(i+1) * fx + Vx(i+1) * fxx
#Quu = luu + fu^T * Vxx(i+1) * fu + Vx(i+1) * fuu
#Qux = lux + fu^T * Vxx(i+1) * fx + Vx(i+1) * fux

class QInfo(object):
    Qx = None
    Qu = None
    Qxx = None
    Quu = None
    Qux = None
    def calc_lx(self, x, u, L):

    def calc(self, x, u, L, f, vinfo):
        lx = calc_lx(x, u, L)
        lu = calc_lu(x, u, L)
        lxx = calc_lxx(x, u, L)
        luu = calc_luu(x, u, L)

class VInfo(object):
    deltaV = None
    Vx = None
    Vxx = None

class DDP(object):
    xgoal = None
    def __init__(self, A, B, Wx, Wu, Wgoal):
        self.A = A
        self.B = B
        self.Wx = Wx
        self.Wu = Wu
        self.Wgoal = Wgoal

    def set_initial_trajectory(self, xlist):
        self.xlist = xlist

    def set_goal(self, xgoal):
        self.xgoal = xgoal

    def costfunctionL(self, x, u):
        return 0.5 * (x.T * Wx * x) + 0.5 * (u.T * Wu * u)

    def costFunctionLgoal(self, x):
        return 0.5 * (x - xgoal).T * Wgoal * (x - xgoal)

    def systemFunctionF(self, x, u):
        return A * x + B * u

    def calcBackwardPath(self):


