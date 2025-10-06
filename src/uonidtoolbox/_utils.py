
import numpy as np
import scipy
import uonidtoolbox as unit
import warnings
import os

import sys
eps = sys.float_info.epsilon


def length(o):
    if isinstance(o, np.ndarray):
        l = np.max(o.shape)
    elif isinstance(o, (int, float, np.int64, np.float64)):
        l = 1
    else:
        l = len(o)
    #endif
    return l
#endfunction

def isempty(obj):
    if isinstance(obj, np.ndarray):
        return True if obj.size == 0 else False
    elif isinstance(obj, (list, dict, unit.struct)):
        return True if len(obj) == 0 else False
    else:
        return False
    #endif
#endfunction


def udisp(msg, guiRunning=0, guiHandle=None):
    if guiRunning:
        #send message to GUI
        print(msg)
    else:
        print(msg)
    #endif
#endfunction

def uwarning(msg):
    warnings.warn(msg)
#endfunction

# def path_to_pkg():
#     # path\to\UNITpy
#     return os.path.abspath(os.path.join(os.path.dirname(unit.__file__), os.pardir, os.pardir))
# #endfunction


def blockhankel(x, nr):
    nx = x.shape[0]
    N = x.shape[1]
    X = np.empty([nx*nr, N-nr+1])
    for idx in range(0, nx):
        X[idx::nx, :] = scipy.linalg.hankel(x[idx, 0:nr], x[idx, nr-1:N])
    #endfor

    return X
#endfunction


# https://math.stackexchange.com/questions/2739271/similarity-transformation-between-two-state-space-ss-models-of-the-same-system
def getSimilarityTransform(SS_1, SS_2):
    nx      = SS_1.A.shape[0]
    ny,nu   = SS_1.D.shape
    I       = np.eye(nx,nx)

    M = np.vstack([
        np.kron(I, SS_2.A) - np.kron(SS_1.A.transpose(), I),
        np.kron(SS_1.B.transpose(), I),
        np.kron(I, SS_2.C)
    ])

    v = np.vstack([
        np.zeros([nx*nx, 1]),
        SS_2.B.reshape(nx*nu, 1),
        SS_1.C.reshape(nx*ny, 1)
    ])

    P = (np.linalg.pinv(M) @ v).reshape(nx,nx)

    # Check residuals are ~= 0
    tol = 1e-16
    r = (v - M @ P.reshape(nx*nx,1)).flatten()
    if r.dot(r) > tol:
        raise Exception(f"residuals of similarity transform were above tolerance: {r.dot(r)} > {tol}")
    #endif

    return P
#enddef


# TODO: implement handling of MIMO, MISO, SIMO systems
from numpy.polynomial.polynomial import polyval as np_polyval
import copy
def m2f(M):
    if 'finishM' not in M: M = unit.startM(M)

    G = copy.deepcopy(M)

    if M.op == 'q':
        ww      = np.exp(  1j*M.w*M.T)
        pdel    = np.exp((-1j*M.w*M.T)*G.delay) # phase lag due to delays on inputs
    else:
        raise Exception("M.op == "+str(M.op)+" not implemented")
    #endif

    match G.type:
        case "ss":
            # TODO: implement handling of SS systems
            raise Exception("m2f() not yet implemented for G.type == ss")
        case _: # polynomial model
            if M.type in ['ar','arx','arma','armax']: G.D = copy.deepcopy(G.A)
            if M.type in ['ar','arma']: G.C = np.array([1.0])
            if M.type in ['ar','arma']: G.B = np.array([0.0])
            if M.type in ['fir']: G.A = np.array([1.0])
            if M.type in ['fir']: G.C = np.array([1.0])
            if M.type in ['fir']: G.D = copy.deepcopy(G.A)

            # TODO: repair dimensions of A,B,C,D polynomials (unit.startM() ?)
            tmpG = unit.struct()
            for p in ['A','B','C','D']:
                if G[p].ndim > 1:
                    tmpG[p] = G[p][0,:]
                else:
                    tmpG[p] = G[p]
                #endif
            #endfor

            A = tmpG.A
            B = tmpG.B

            G.G = np_polyval(1/ww, B)/np_polyval(1/ww, A)
            pp = pdel
            G.G = G.G * pp # pp is responsible for the effect of delay on input

            C = np.hstack([np.zeros(tmpG.D.size - tmpG.C.size), tmpG.C])
            D = tmpG.D
            G.H = np_polyval(1/ww, C)/np_polyval(1/ww, D)
    #endmatch

    # For time series case, noise spectral factor masquerades as dynamic freq resp.
    if G.nu<1: G.G = G.H

    return G
#endfunction



def theta2m(theta, M):
    idx = 0
    M.B = theta[idx:idx+M.nB[0]+1]; idx += M.nB[0]+1
    M.A = theta[idx:idx+M.nA[0]]; M.A = np.hstack([1, M.A]); idx += M.nA[0]

    return M
#endfunction


def m2theta(M):
    # TODO: Use information about model structure to determine which objects should be loaded into theta vector
    # Can be either polynomials or state space matrices

    # In case of SISO OE model type, load B(q) and A(q) polynomials into theta
    # First, determine if initial polynomial estimates are in M, otherwise they are orders
    Mt = unit.struct()
    for p in ['B', 'A']:
        if isinstance(M[p], np.ndarray):
            if M[p].size > 1: # its a polynomial
                Mt[p] = M[p]
            else: # its an order
                Mt[p] = np.ones(M[p].ravel()[0])
            #endif
        elif isinstance(M[p], (int, float)):
            Mt[p] = np.ones(M[p].ravel()[0])
        else:
            raise Exception(f"M.{p} is neither int nor np.ndarray")
        #endif
    #endfor

    M.nA = np.array([Mt.A.size - 1])
    M.nB = np.array([Mt.B.size - 1])

    # Stack polynomials in vector
    theta = np.hstack([Mt.B.ravel(), Mt.A.ravel()[1:]])

    return theta
#endfunction
