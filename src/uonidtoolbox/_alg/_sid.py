
import uonidtoolbox as unit
import numpy as np
import scipy
import copy


def sid(Z, M=unit.struct(), OPT=unit.struct()):

    # Extract inputs and outputs specified
    y,u,ny,nu,N = unit._setup._startZ._Z2data(Z)

    # Check which parts of model structure were unspecified and set to defaults.
    if not M:
        unit._utils.uwarning('Need to specify initial model structure M!')
    #endif
    
    M.type = 'ss'
    M = unit._setup.startM(Z,M)
    order = M.nx

    # Start with G = M
    G = copy.deepcopy(M)

    # Include delays specified in model structure on inputs
    for r in range(0, nu):
        u[:,r] = np.hstack([np.zeros(M.delay[r]), u[0:N-M.delay[r], r]])
    #endfor

    # Check what options not specified explicitly by user and then set to defaults
    OPT = unit._setup.startOPT(OPT)
    if 'bw' not in OPT: OPT.bw = 0.5/M.T
    if 'horizon' not in OPT:
        OPT.horizon = np.min(np.floor( np.array([2.5*order, (N+1-ny)/(2*nu+ny+2)]) ))
    #endif
    if 'alg' not in OPT: 
        OPT.alg = 'n4sid'
    elif OPT.alg.lower() == 'sid':
        OPT.alg = 'n4sid'
    elif OPT.alg.lower() not in ['n4sid', 'cca']:
        OPT.alg = 'n4sid'
    #endif

    # # Forward and backward horizons are the same
    # f = OPT.horizon
    # p = f
    i = OPT.horizon # = p = f

    # Form block Hankel Matrix from inputs and outputs
    u_hankel = unit._utils.blockhankel(u.transpose(), 2*i)
    y_hankel = unit._utils.blockhankel(y.transpose(), 2*i)

    Up = u_hankel[:nu*i , :] # U_{0: i-1}
    Uf = u_hankel[ nu*i:, :] # U_{i:2i-1}
    Yp = y_hankel[:ny*i , :] # Y_{0: i-1}
    Yf = y_hankel[ ny*i:, :] # Y_{i:2i-1}

    H = np.vstack([Up, Uf, Yp, Yf])

    R = np.linalg.qr(H.transpose(), mode='r').transpose() # lower triangular qr factorisation
    R = np.tril(R)

    R_22_14 = _extract_block_from_R(R, 2, 2, 1, 4, nu, ny, i)
    R_55_14 = _extract_block_from_R(R, 5, 5, 1, 4, nu, ny, i)

    R_56_14 = _extract_block_from_R(R, 5, 6, 1, 4, nu, ny, i)
    R_66_15 = _extract_block_from_R(R, 6, 6, 1, 5, nu, ny, i)

    R_14_14 = _extract_block_from_R(R, 1, 4, 1, 4, nu, ny, i)
    R_15_15 = _extract_block_from_R(R, 1, 5, 1, 5, nu, ny, i)
    R_15_14 = _extract_block_from_R(R, 1, 5, 1, 4, nu, ny, i)

    term_i      = R_56_14 @ np.linalg.pinv(R_14_14)
    term_ip1    = R_66_15 @ np.linalg.pinv(R_15_15)

    L1_0_L3_i   = np.hstack([term_i  [:,0:nu*i],   np.zeros([ny*i  ,nu*i]),   term_i  [:,2*nu*i:]])
    L1_0_L3_ip1 = np.hstack([term_ip1[:,0:nu*i+1], np.zeros([ny*i-1,nu*i-1]), term_ip1[:,2*nu*i:]])

    O, Sigma, V = np.linalg.svd(L1_0_L3_i @ R_14_14)

    O_d = O[:, :order]
    s_d = Sigma[:order]

    X_i     = np.diag(np.sqrt(1/s_d)) @ O_d.transpose() @ L1_0_L3_i @ R_14_14
    X_ip1   = np.diag(np.sqrt(1/s_d)) @ np.linalg.pinv(O_d[:-ny, :]) @ L1_0_L3_ip1 @ R_15_14

    # \mathcal{L} = [X_ip1; Y_ii] @ pinv([X_i; U_ii])
    ABCD = np.vstack([X_ip1, R_55_14]) @ np.linalg.pinv(np.vstack([X_i, R_22_14]))

    G.ss = unit.struct()
    G.ss.A = ABCD[0:order, 0:order].reshape(order, order)
    G.ss.B = ABCD[0:order, order:order+nu].reshape(order, nu)
    G.ss.C = ABCD[order:order+ny, 0:order].reshape(ny, order)
    G.ss.D = ABCD[order:order+ny, order:order+nu].reshape(ny, nu)

    # residuals
    rho = np.vstack([X_ip1, R_55_14]) - ABCD @ np.vstack([X_i, R_22_14])

    rho1 = rho[0:order]           # states
    rho2 = rho[order:order+ny]    # outputs

    G.ss.Q = ((1/N)*(rho1 @ rho1.transpose())).reshape(order, order)
    G.ss.R = ((1/N)*(rho2 @ rho2.transpose())).reshape(ny, ny)
    G.ss.S = ((1/N)*(rho1 @ rho2.transpose())).reshape(order, ny)

    G.nx = order
    G.th = ABCD


    return G
#endfunction


# R_{5:6, 1:4} = _extract_block_from_R(R, 5, 6, 1, 4)
def _extract_block_from_R(R, brs, bre, bcs, bce, nu, ny, i):
    block_dims = np.array([ nu*i, nu, nu*(i-1), ny*i, ny, ny*(i-1) ])

    rs = np.int64(block_dims.dot(np.hstack([np.ones(brs-1), np.zeros(6-brs+1)])))
    re = np.int64(block_dims.dot(np.hstack([np.ones(bre  ), np.zeros(6-bre  )])))
    cs = np.int64(block_dims.dot(np.hstack([np.ones(bcs-1), np.zeros(6-bcs+1)])))
    ce = np.int64(block_dims.dot(np.hstack([np.ones(bce  ), np.zeros(6-bce  )])))

    return R[rs:re, cs:ce]
#endfunction

