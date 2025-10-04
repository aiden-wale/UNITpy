
import uonidtoolbox as unit
import numpy as np


def rosenbrock(x, y, compute_gradient=False):
    a = 1; b = 100

    amx     = a - x
    ymx2    = y - x*x

    f = amx*amx + b*ymx2*ymx2

    if not compute_gradient:
        return f
    #endif

    g = np.array([
            -2*amx - 4*b*ymx2,
            2*b*ymx2
        ])

    return f,g
#endfunction


def expcurve(th, t, compute_gradient=False):
    f = (th[0] * np.exp(th[1]*t) + th[2]).reshape(t.size,1)

    if not compute_gradient:
        return f
    #endif

    # df/dtheta
    J = np.ndarray([t.size, th.size])
    J[:,0] = np.exp(th[1]*t)
    J[:,1] = th[0] * t * np.exp(th[1]*t)
    J[:,2] = 1.0
    return f,J
#enddef


