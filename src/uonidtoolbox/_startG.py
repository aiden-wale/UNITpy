
import uonidtoolbox as unit
import numpy as np
import scipy
import copy


def startG(Z,M,OPT):
    OPT = unit.startOPT(OPT)
    Z   = unit.startZ(Z)

    if M.type in ['oe', 'bj']:
        g_SM = _SteiglitzMcBride_Init(Z, M, OPT)
    #endif

    M.A = g_SM.A
    M.B = g_SM.B

    return M
#endfunction


def _SteiglitzMcBride_Init(Z, M, OPT):
    if OPT.dsp: unit._utils.udisp("Finding initial dynamics model via Steiglitz-McBride...")

    ZZ = copy.deepcopy(Z)
    MM = copy.deepcopy(M)
    # MM.delay[:] = 0

    # Create intial guess via arx
    g = unit._alg.barx(ZZ,MM,OPT)
    # TODO: Make sure the A polynomial is stable

    gfav = copy.deepcopy(g)
    itfav = 0

    th      = unit._alg._gn._m2theta(g)
    cost    = unit._alg._gn._VN(th, Z, MM, OPT, compute_gradient=False)
    if OPT.dsp: unit._utils.udisp(f"iter#: {0:<3} |  cost = {cost:<10.5e}")

    for it in range(1, OPT.smits+1):
        # Perform prefiltering step
        uh = scipy.signal.lfilter(1, g.A, Z.u.transpose()).transpose()
        yh = scipy.signal.lfilter(1, g.A, Z.y.transpose()).transpose()

        # Solve for polynomials using prefiltered data
        ZZ.u = uh.copy()
        ZZ.y = yh.copy()
        g = unit._alg.barx(ZZ, MM, OPT)
        # TODO: Make sure the A polynomial is stable

        th      = unit._alg._gn._m2theta(g)
        costnew = unit._alg._gn._VN(th, Z, MM, OPT, compute_gradient=False)
        if OPT.dsp: unit._utils.udisp(f"iter#: {it:<3} |  cost = {costnew:<10.5e}")

        if costnew<cost:
            cost    = costnew
            gfav    = g
            itfav   = it
        #endif
    #endfor
    g = gfav
    if OPT.dsp: unit._utils.udisp(f"Best results achieved at iteration {itfav} with cost = {cost:<10.5e}.\n")

    return g
#endfunction