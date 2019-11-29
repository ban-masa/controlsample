import numpy as np

#Qx = lx + fx^T * nextVx(i+1)
#Qu = lu + fu^T * nextVx(i+1)
#Qxx = lxx + fx^T * Vxx(i+1) * fx + Vx(i+1) * fxx
#Quu = luu + fu^T * Vxx(i+1) * fu + Vx(i+1) * fuu
#Qux = lux + fu^T * Vxx(i+1) * fx + Vx(i+1) * fux

#deltaV = -0.5 * k^T * Quu * k
#Vx = Qx - K^T * Quu * k
#Vxx = Qxx - K^T * Quu * K

class QInfo(object):
    Qx = None
    Qu = None
    Qxx = None
    Quu = None
    Qux = None
    x = None
    xref = None
    u = None
    L = None
    def __init__(self):
        return

    def set_params(self, x_, xref_, u_, L_):


    def calc_lx(self, x, xref, u, L):
        dim = x.shape[0]
        ret = np.zeros((dim, 1))
        for i in range(dim):
            xtemp = x.copy()
            xtemp[i][0] = xtemp[i][0] + 0.00001
            ret[i][0] = (L(xtemp, u) - L(x, u)) / 0.00001
        return ret.T[0]

    def calc_lu(self, x, u, L):
        dim = u.shape[0]
        ret = np.zeros((dim, 1))
        for i in range(dim):
            utemp = u.copy()
            utemp[i][0] = utemp[i][0] + 0.00001
            ret[i][0] = (L(x, utemp) - L(x, u)) / 0.00001
        return ret.T[0]

    def calc_lxx(self, x, u, L):
        dim = x.shape[0]
        ret = np.zeros((dim, dim))
        for i in range(dim):
            xtemp = x.copy()
            xtemp[i][0] = xtemp[i][0] + 0.00001
            ret[:,i] = (calc_lx(xtemp, u, L) - calc_lx(x, u, L)) / 0.00001
        return ret

    def calc_lxu(self, x, u, L):
        xdim = x.shape[0]
        udim = u.shape[0]
        ret = np.zeros((xdim, udim))
        for i in range(udim):
            utemp = u.copy()
            utemp[i][0] = utemp[i][0] + 0.00001
            ret[:,i] = (calc_lx(x, utemp, L) - calc_lx(x, u, L)) / 0.00001
        return ret

    def calc_lux(self, x, u, L):
        xdim = x.shape[0]
        udim = u.shape[0]
        ret = np.zeros((udim, xdim))
        for i in range(xdim):
            xtemp = x.copy()
            xtemp[i][0] = xtemp[i][0] + 0.00001
            ret[:,i] = (calc_lu(xtemp, u, L) - calc_lu(x, u, L)) / 0.00001
        return ret

    def calc_luu(self, x, u, L):
        dim = u.shape[0]
        ret = np.zeros((dim, dim))
        for i in range(dim):
            utemp = u.copy()
            utemp[i][0] = utemp[i][0] + 0.00001
            ret[:,i] = (calc_lu(x, utemp, L) - calc_lu(x, u, L)) / 0.00001
        return ret

    def calc_fx(self, x, u, f):
        dim = x.shape[0]
        ret = np.zeros((dim, dim))
        for i in range(dim):
            xtemp = x.copy()
            xtemp[i][0] = xtemp[i][0] + 0.00001
            ret[:,i] = (f(xtemp, u).T - f(x, u).T) / 0.00001
        return ret
    
    def calc_fu(self, x, u, f):
        xdim = x.shape[0]
        udim = u.shape[0]
        ret = np.zeros((xdim, udim))
        for i in range(udim):
            utemp = u.copy()
            utemp[i][0] = utemp[i][0] + 0.00001
            ret[:,i] = (f(x, utemp).T - f(x, u).T) / 0.00001
        return ret

    def calc(self, x, xref, u, L, f, vinfo):
#Qx = lx + fx^T * nextVx(i+1)
#Qu = lu + fu^T * nextVx(i+1)
#Qxx = lxx + fx^T * Vxx(i+1) * fx + Vx(i+1) * fxx
#Quu = luu + fu^T * Vxx(i+1) * fu + Vx(i+1) * fuu
#Qux = lux + fu^T * Vxx(i+1) * fx + Vx(i+1) * fux
        lx = calc_lx(x, u, L).T
        lu = calc_lu(x, u, L).T
        lxx = calc_lxx(x, u, L)
        luu = calc_luu(x, u, L)
        lux = calc_lux(x, u, L)
        lxu = calc_lxu(x, u, L)
        fx = calc_fx(x, u, f)
        fu = calc_fu(x, u, f)
        Qx = lx + np.dot(fx.T, vinfo.Vx)
        Qu = lu + np.dot(fu.T, vinfo.Vx)
        Qxx = lxx + np.dot(np.dot(fx.T, vinfo.Vxx), fx)
        Quu = luu + np.dot(np.dot(fu.T, vinfo.Vxx), fu)
        Qux = lux + np.dot(np.dot(fu.T, vinfo.Vxx), fx)

    def update_u(self, 

class VInfo(object):
#deltaV = -0.5 * k^T * Quu * k
#Vx = Qx - K^T * Quu * k
#Vxx = Qxx - K^T * Quu * K
    V = None
    deltaV = None
    Vx = None
    Vxx = None
    def __init__(self):
        return

    def init(self):
        return

    def update(self):
        return


class DDP(object):
    xgoal = None
    QInfo qinfo = None
    VInfo vinfo = None
    xlist = None
    ulist = None
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

    def costfunctionL(self, x, xref, u):
        return 0.5 * (x.T * Wx * x) + 0.5 * (u.T * Wu * u)

    def costFunctionLgoal(self, x):
        return 0.5 * (x - xgoal).T * Wgoal * (x - xgoal)

    def systemFunctionF(self, x, u):
        return A * x + B * u

#初期値はどうやって作るのだろうか。初期入力列から初期状態列を作る？
    def calcBackwardPath(self):
        vinfo.init()
        for xu in self.xulist_revert:
            qinfo.calc(x, u, L, f, vinfo)

            vinfo.calc(qinfo)

        
        return

    def calcForwardPath(self):

        return


