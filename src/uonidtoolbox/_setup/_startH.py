
import uonidtoolbox as unit
import numpy as np
import scipy
import copy


def startH(Z,M,OPT):
    OPT = unit._setup.startOPT(OPT)
    Z   = unit._setup.startZ(Z)

    if M.type in ['oe']:
        M.C  = np.array([1])
        M.D  = np.array([1])
    elif M.type in ['bj']:
        g_HR = _HannanRissanen_Init(Z, M, OPT)
        M.C  = g_HR.C
        M.D  = g_HR.D
    #endif

    return M
#endfunction


def _HannanRissanen_Init(Z, M, OPT):
    if OPT.dsp: unit._utils.udisp("Finding initial noise model via Hannan-Rissanen...")
    
    y,u,ny,nu,Ny = unit._setup._startZ._Z2data(Z)

    # Evaluate dynamics term w = G*u, using existing estimate
    w = scipy.signal.lfilter(M.B, M.A, u.transpose())

    # Fit high order AR model to coloured noise term: D*v = (y - G*u)
    ZZ = copy.deepcopy(Z)
    MM = copy.deepcopy(M)

    v = y.transpose() - w
    ZZ.nu = 0
    ZZ.y  = v.reshape(v.size,1)
    ZZ.Ny = y.size
    MM.nA = np.min([100, Ny/10])
    MM.A  = MM.nA
    MM.nB = 1
    MM.B  = np.array([1])
    g = unit._alg.barx(ZZ,MM,OPT)

    # Estimate e by filtering through D. Ignore transients from first <order> data points
    e = scipy.signal.lfilter(g.A, 1, v)
    e = e[:, g.A.size-1:]
    v = v[:, g.A.size-1:]

    # Fit ARX model to v = C/D*e
    ZZ = copy.deepcopy(Z)
    MM = copy.deepcopy(M)

    ZZ.u  = e.reshape(e.size,1)
    ZZ.y  = v.reshape(v.size,1)
    ZZ.Ny = v.size
    MM.nB = MM.nC
    MM.nA = MM.nD
    MM.B  = MM.C
    MM.A  = MM.D
    g = unit._alg.barx(ZZ,MM,OPT)

    M.C = g.B
    M.D = g.A

    # gfav = copy.deepcopy(g)
    # itfav = 0

    # th      = unit._utils.m2theta(g)
    # cost    = unit._optimisation.VN(th, Z, MM, OPT, compute_gradient=False)
    # if OPT.dsp: unit._utils.udisp(f"iter#: {0:<3} |  cost = {cost:<10.5e}")

    # for it in range(1, OPT.smits+1):
    #     # Perform prefiltering step
    #     uh = scipy.signal.lfilter(1, g.A, Z.u.transpose()).transpose()
    #     yh = scipy.signal.lfilter(1, g.A, Z.y.transpose()).transpose()

    #     # Solve for polynomials using prefiltered data
    #     ZZ.u = uh.copy()
    #     ZZ.y = yh.copy()
    #     g = unit._alg.barx(ZZ, MM, OPT)
    #     # TODO: Make sure the A polynomial is stable

    #     th      = unit._utils.m2theta(g)
    #     costnew = unit._optimisation.VN(th, Z, MM, OPT, compute_gradient=False)
    #     if OPT.dsp: unit._utils.udisp(f"iter#: {it:<3} |  cost = {costnew:<10.5e}")

    #     if costnew<cost:
    #         cost    = costnew
    #         gfav    = g
    #         itfav   = it
    #     #endif
    # #endfor
    # g = gfav
    # if OPT.dsp: unit._utils.udisp(f"Best results achieved at iteration {itfav} with cost = {cost:<10.5e}.\n")

    if OPT.dsp: unit._utils.udisp(f"Best results achieved at iteration {0} with cost = {0.0:<10.5e}.\n")

    return M
#endfunction