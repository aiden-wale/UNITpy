
import uonidtoolbox as unit
import numpy as np
import scipy
import copy


def barx(Z, M=unit.struct(), OPT=unit.struct()):
    Z = unit.startZ(Z)
    y,u,ny,nu,Ny,Z = unit._startZ._Z2data(Z)

    # Unspecified parts of OPT -> defaults
    OPT = unit.startOPT(OPT)
    if 'alg' not in OPT:
        OPT.alg = unit.struct()
        if 'type' not in OPT.alg:
            OPT.alg.type = 'block'
        #endif
    #endif

    # Unspecified parts of M -> defaults
    if OPT.n >= Ny:
        raise Exception("Cannot have OPT.n larger than height of Z!")
    #endif

    if unit._utils.isempty(M):
        raise Exception("Need to specify initial model structure M!")
    else:
        M = unit.startM(Z,M)
    #endif

    # Include delays specified in model structure on inputs
    for r in range(0,nu):
        u[:,r:r+1] = np.vstack([ np.zeros([M.delay[r], 1]) , u[0:Ny-M.delay[r], r:r+1] ])
    #endfor

    if nu>0: # ARX
        m = np.sum(M.nB) + nu
    else: # AR
        m = 0
    #endif

    n = np.max(M.nA)

    # # TODO: Pass input through a non-linearity if required by model structure
    # x = unit.u2x(u,M)
    x = u

    # Form regressor matrix for input
    PHI = np.empty([Ny, n+m])

    # Insert u data into PHI
    idx = 0
    for r in range(0,nu):
        # PHIu = toeplitz(u, [u(1), zeros(1,M.nB)])
        PHI[:, idx:idx+M.nB[r]+1] = scipy.linalg.toeplitz(x[:,r], np.hstack([x[0,r], np.zeros(M.nB[r])]))
        idx += M.nB[r]+1
    #endif

    # Insert y data into PHI
    PHI[0,  m:m+n] = np.zeros(n)
    PHI[1:, m:m+n] = scipy.linalg.toeplitz(-y[0:-1,0], np.hstack([-y[0,0], np.zeros(n-1)]))



    # Save initial model into G
    G = copy.deepcopy(M)

    # Get model estimate via least squares
    G.th = np.linalg.lstsq(PHI[OPT.n:Ny, :], y[OPT.n:Ny])[0]

    # Unpack G.th into A,B
    # TODO: handle MISO
    G.B = G.th[0:m, 0]
    G.A = np.hstack([1, G.th[m:m+n,0]])

    G.phi = PHI
    
    return G
#endfunction

