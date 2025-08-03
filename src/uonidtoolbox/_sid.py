
import numpy as np
# import scipy
import uonidtoolbox as unit
import copy


def sid(Z,M={},OPT={}):

    # TODO: complete this function
    raise Exception("sid: not implemented yet")
    return 0

    # Extract inputs and outputs specified
    y,u,ny,nu,N,Z = unit.Z2data(Z)

    # Check which parts of model structure were unspecified and set to defaults.
    if not M:
        unit._utils.uwarning('Need to specify initial model structure M!')
    #endif
    
    M['type'] = 'ss'
    M = unit.startM(Z,M)
    order = M['nx']

    # Start with G = M
    G = copy.deepcopy(M)

    # Include delays specified in model structure on inputs
    for r in range(0, nu):
        u[:,r] = np.vstack([np.zeros([M['delay'][r], 1]), u[0:N-M['delay'][r], r]])
    #endfor

    # Check what options not specified explicitly by user and then set to defaults
    OPT = startOPT(OPT)
    if 'bw' not in OPT: OPT['bw'] = 0.5/M['T']
    if 'horizon' not in OPT:
        OPT['horizon'] = np.min(np.floor( np.array([2.5*order, (N-OPT['n']+1-ny)/(2*nu+ny+2)]) ))
    #endif
    if 'alg' not in OPT: 
        OPT['alg'] = 'n4sid'
    elif OPT['alg'].lower() == 'sid':
        OPT['alg'] = 'n4sid'
    elif OPT['alg'].lower() not in ['n4sid', 'cca']:
        OPT['alg'] = 'n4sid'
    #endif

    # Forward and backward horizons are the same
    f = OPT['horizon']
    p = f

    # Form block Hankel Matrices from inputs and outputs
    Y = np.empty([     ny*f, N-(p+f)-OPT['n']+1])
    U = np.empty([     nu*f, N-(p+f)-OPT['n']+1])
    Z = np.empty([(nu+ny)*p, N-(p+f)-OPT['n']+1])

    if M['op'] == 'q':
        for k in range(0, f):
            Y[k*ny:(k+1)*ny,:] = y[OPT['n']+p+k:N-f+k+1,0:ny].T
            if nu>0:
                U[(f-k-1)*nu:(f-k)*nu,:] = u[OPT['n']+p+k:N-f+k+1,0:nu].T
            #endif
        #endfor
        for k in range(0, p):
            Z[nu*p+k*ny:nu*p+(k+1)*ny,:] = y[OPT['n']+p-k-1:N-f-k,0:ny].T
            if nu>0:
                Z[k*nu:(k+1)*nu,:] = u[OPT['n']+p-k-1:N-f-k,0:nu].T
            #endif
        #endfor
    elif M['op'] == 'd':
        # TODO: implement 'd' operator
        unit._utils.uwarning("sid: M['op'] == 'd': 'd' operator is not implemented yet")
    #endif

    
    # Switch to ortho basis for quantities involved - much more efficient flop wise
    R = np.linalg.qr(np.hstack([U.T, Z.T, Y.T]), mode='r')
    R = np.triu(R)
    R11     = R[0:f*nu, 0:f*nu]
    R12     = R[0:f*nu, f*nu:(f+p)*nu+p*ny]
    R12p    = R[0:(f-1)*nu, (f-1)*nu:(f+p)*nu+(p+1)*ny]
    R13     = R[0:f*nu, (f+p)*nu+p*ny:(f+p)*(nu+ny)]
    R22     = R[f*nu:(f+p)*nu+p*ny, f*nu:(f+p)*nu+p*ny]
    R22p    = R[(f-1)*nu:(f+p)*nu+(p+1)*ny, (f-1)*nu:(f+p)*nu+(p+1)*ny]
    R23     = R[f*nu:(f+p)*nu+p*ny, (f+p)*nu+p*ny:(f+p)*(nu+ny)]
    R23p    = R[(f-1)*nu:(f+p)*nu+(p+1)*ny, (f+p)*nu+(p+1)*ny:(f+p)*(nu+ny)]
    R33     = R[(f+p)*nu+p*ny:, (f+p)*nu+p*ny:(f+p)*(nu+ny)]
    # R33p    = R[(f+p)*nu+(p+1)*ny:, (f+p)*nu+(p+1)*ny:(f+p)*(nu+ny)]

    # Linear regression solution w.r.t to orthonormal bases
    beta    = R23.T*R22*np.linalg.pinv(R22.T*R22)
    bplus   = R23p.T*R22p*np.linalg.pinv(R22p.T*R22p)
    Zp      = np.hstack([R12.T, R22.T])  # Y*PI*Z'*(Z*PI*Z')^(-1) = beta
    Zplus   = np.hstack([R12p.T, R22p.T]) # For Xplus computation

    # Compute weightings according to algorithm variant
    if OPT['alg'].lower() == 'cca':
        U,S,V = np.linalg.svd

    # Estimate subspace spanned by columns of extended observabililty matrix O

    # Extract estimated states - used for Q matrix estimation + par est in n4sid

    # Extract estimates of A,B,C,D from the estimate of the observability matrix O and the estimate of the regressor beta.

    # Compute state and noise covariances from residuals

    # Convert from state space to transfer function form firstly remove G.A and G.B in case they came with the initial model

    # Return final estimate in innovations form - suppress warnings/errors from dare

    # Convert from state space to transfer function form

    # Fill in last remaining details about the estimated model structure

    # Fill in bilinear info

    # Add legend for prospective plotting

    # Add in estimation of innovations variance

    # Record that VNss should be used to compute prediction erros for validation

    # Finally, make sure that M.theta reflects the model



    return G
#endfunction

