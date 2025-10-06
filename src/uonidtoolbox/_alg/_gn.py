
import uonidtoolbox as unit
import numpy as np
import scipy
import copy


def gn(Z, M=unit.struct(), OPT=unit.struct()):

    # Extract inputs and outputs specified
    y,u,ny,nu,N,Z = unit._startZ._Z2data(Z)

    M = unit.startM(M)

    # for k in ['D', 'K', 'F', 'G', 'X1']:
    #     if 'est'+k in M:
    #         if not M['est'+k]:
    #             M.ss[k] = np.array([[]])
    #         #endif
    #     #endif
    # #endif

    # Construct parameter vector from system polynomials
    theta0 = _m2theta(M)

    G = copy.deepcopy(M)

    # Create handle to cost function
    costfcn = lambda th, compute_gradient=False: _VN(th, Z, M, OPT, compute_gradient=compute_gradient)

    # Call the gn routine
    thetastar = _gn_impl(costfcn, theta0, maxiter=500, disp=OPT.dsp)

    G.th = thetastar
    G = _theta2m(G.th, G)

    return G
#endfunction


def _gn_impl(costfcn, theta0, maxiter=100, disp=1):
    c1 = 1e-4   # Wolf constant c1
    c2 = 0.9    # Wolf constant c2

    if disp: unit._utils.udisp("Beginning system estimation via Gauss-Newton...")
    f = costfcn(theta0, compute_gradient=False)
    if disp: unit._utils.udisp(f"iter#: {0:<4} |  cost = {f:<10.5e} |  Newton decrement = {'':>12} |  alpha = ")

    if not np.all(np.isfinite(theta0)): raise Exception("Initial parameter estimate is not finite.")

    # Do iterations until solution is found or maxiter is reached
    theta = theta0.copy()
    for it in range(1, maxiter+1):
        # Get cost and gradient of objective function
        f,pe,g,J = costfcn(theta, compute_gradient=True)

        if not np.isfinite(f): raise Exception("Cost is not finite at initial parameter estimate.")

        p   = _computeMinimizingDirection(J,pe)
        pTg = p.ravel().dot(g.ravel())

        # Check stopping criteria
        if np.abs(pTg) < 1e-12:
            if disp: unit._utils.udisp(f"Local minimum found after {it} iteration(s) with cost: {f:<10.5e}\n")
            break
        #endif

        # Find a step length which satisfies Wolfe condition
        alpha = 1.0
        for i in range(0, 51):
            f_n = costfcn(theta + alpha*p.ravel(), compute_gradient=False)
            if f_n < f + c1*alpha*pTg: # Wolfe condition 1
                break
                # f_n,_,g_n,_ = costfcn(theta + alpha*p.ravel(), compute_gradient=True)
                # pTg_n = p.ravel().dot(g_n.ravel())
                # if np.abs(pTg_n) > c2*np.abs(pTg): # Wolfe condition 2
                #     break
                # #endif
            #endif
            alpha /= 2
        #endfor

        # Update theta by suitable step and continue
        theta = theta + alpha*p.ravel()
        if disp: unit._utils.udisp(f"iter#: {it:<4} |  cost = {f_n:<10.5e} |  Newton decrement = {pTg:<10.5e} |  alpha = {alpha:<10.5e}")
    #endfor

    return theta
#endfunction


def _computeMinimizingDirection(J,pe):
    # Compute minimizing direction
    U,s,V = np.linalg.svd(J, full_matrices=False)
    rankJ = len(np.where(s > 1e-14)[0])
    U1 = U[:,0:rankJ]
    s1 = s[0:rankJ]
    V1 = V[0:rankJ,:].transpose()
    p = -V1 @ ((U1.transpose() @ pe) / s1.reshape(s1.size, 1))

    return p
#endfunction


def _theta2m(theta, M):
    idx = 0
    M.B = theta[idx:idx+M.nB[0]+1]; idx += M.nB[0]+1
    M.A = theta[idx:idx+M.nA[0]]; M.A = np.hstack([1, M.A]); idx += M.nA[0]

    return M
#endfunction


def _m2theta(M):
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


def _VN(theta, Z, M, OPT, compute_gradient=False):
    # TODO: handle systems other than OE
    # TODO: handle MISO polynomial case(s)
    # Extract inputs and outputs specified
    y,u,ny,nu,Ny,Z = unit._startZ._Z2data(Z)

    # Include delays specified in model structure on inputs
    for r in range(0,nu):
        u[:,r:r+1] = np.vstack([ np.zeros([M.delay[r], 1]) , u[0:Ny-M.delay[r], r:r+1] ])
    #endfor

    # Get polynomials in model structure form, from theta
    Mn = _theta2m(theta, M)
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


