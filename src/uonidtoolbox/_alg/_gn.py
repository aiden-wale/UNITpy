
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
    theta0 = unit._utils.m2theta(M)

    G = copy.deepcopy(M)

    # Create handle to cost function
    costfcn = lambda th, compute_gradient=False: unit._objective.VN(th, Z, M, OPT, compute_gradient=compute_gradient)

    # Call the gn routine
    thetastar = _gn_impl(costfcn, theta0, maxiter=500, disp=OPT.dsp)

    G.th = thetastar
    G = unit._utils.theta2m(G.th, G)

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

        # Compute minimizing direction
        U,s,V = np.linalg.svd(J, full_matrices=False)
        rankJ = len(np.where(s > 1e-14)[0])
        U1 = U[:,0:rankJ]
        s1 = s[0:rankJ]
        V1 = V[0:rankJ,:].transpose()
        p = -V1 @ ((U1.transpose() @ pe) / s1.reshape(s1.size, 1))

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

