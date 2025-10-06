
import uonidtoolbox as unit
import numpy as np
import scipy
import copy


def VN(theta, Z, M, OPT, compute_gradient=False):
    # TODO: handle systems other than OE
    # TODO: handle MISO polynomial case(s)
    # Extract inputs and outputs specified
    y,u,ny,nu,Ny,Z = unit._startZ._Z2data(Z)

    # Include delays specified in model structure on inputs
    for r in range(0,nu):
        u[:,r:r+1] = np.vstack([ np.zeros([M.delay[r], 1]) , u[0:Ny-M.delay[r], r:r+1] ])
    #endfor

    # Get polynomials in model structure form, from theta
    Mn = unit._utils.theta2m(theta, M)
    a = Mn.A.ravel()
    b = Mn.B.ravel()
    c = Mn.C.ravel()
    d = Mn.D.ravel()

    yh = scipy.signal.lfilter(b, a, u.transpose()).transpose()

    # Compute prediction errors and cost
    pe      = yh - y
    cost    = 0.5*pe.ravel().dot(pe.ravel())/Ny

    if not compute_gradient:
        return cost
    #endif

    # Compute jacobian
    J = np.ndarray([Ny, theta.size])
    idx = 0
    ac = np.convolve(a, c)
    for k in range(0, M.nB[0]+1):
        num = np.zeros(M.nB[0]+1); num[k] = 1
        J[:, idx] = scipy.signal.lfilter(np.convolve(num,d), ac, u.transpose()).ravel()
        idx += 1
    #endfor
    for k in range(1, M.nA[0]+1):
        num = np.zeros(M.nA[0]+1); num[k] = 1
        J[:, idx] = -scipy.signal.lfilter(np.convolve(d,num), ac, yh.transpose()).ravel()
        idx += 1
    #endfor

    g = J.transpose() @ pe/Ny

    return cost,pe,g,J
#enddef

