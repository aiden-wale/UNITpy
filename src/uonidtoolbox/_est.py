
import numpy as np
import uonidtoolbox as unit


def est(Z, M, OPT):

    G = 0

    if not Z:
        raise Exception("Need to specify data (Z)!")
    elif not M:
        Z = unit.startZ(Z)
        m.nx = np.min([20, np.ceil(Z.Ny/10)])
        gsid = unit.subspace(Z, m)
        lsin = length(gsid.sing)
        vv = np.vstack([np.linspace(0, (lsin-1)/lsin, lsin), gsid.sing.T/gsid.sing[0]])
        mi = np.argmin(np.sum(vv*vv, 0))
        mv = np.sum(vv*vv, 0)[mi]
        M.A = mi + 1
        OPT = {}
    elif not OPT:
        OPT = {}
    #endif

    Z   = unit.startZ(Z)
    M   = unit.startM(Z,M)
    OPT = unit.startOPT(OPT,M)
    ep  = unit.estmap(Z,M,OPT)

    if OPT.dsp:
        dblines = "===================================================================="
        unit._utils.udisp("\n" + dblines)
        unit._utils.udisp("START ESTIMATION PROCESS:")
        unit._utils.udisp("Estimating parameters for '" + M.type.upper() + "' model structure using '" + M.op.lower() + "' operator.")
        unit._utils.udisp(ep.modelEquations)
        unit._utils.udisp("INITIALISATION:")
    #endif

    # Init nonlinear parts of model if necessary
    # if not isempty(ep.startNL):
    #     M = 

    # Init estimate of system dynamics if necessary

    # Init noise model if necessary

    
    # Now call appropriate estimation algorithm

    # Fill in components of returned (estimated) model structure


    if OPT.dsp:
        unit._utils.udisp("END ESTIMATION PROCESS")
        unit._utils.udisp(dblines + "\n")
    #endif

    return G













