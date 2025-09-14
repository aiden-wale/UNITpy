
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
    elif isinstance(o, int):
        l = 1
    elif isinstance(o, float):
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
        np.kron(I, SS_2.A) - np.kron(SS_1.A.T, I),
        np.kron(SS_1.B.T, I),
        np.kron(I, SS_2.C)
    ])

    v = np.vstack([
        np.zeros([nx*nx, 1]),
        SS_2.B.reshape(nx*nu, 1),
        SS_1.C.reshape(nx*ny, 1)
    ])

    T = (np.linalg.pinv(M) @ v).reshape(nx,nx)
    return T
#enddef


