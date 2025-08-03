
import numpy as np
import scipy
import uonidtoolbox as unit
import copy


def fir(Z,M={},OPT={}):
    G = 0
    
    Z = unit.startZ(Z)
    y,u,ny,nu,Ny,Z = unit._startZ.Z2data(Z)

    # Unspecified parts of OPT -> defaults
    OPT = unit.startOPT(OPT)
    if 'type' not in OPT['alg']:
        OPT['alg'] = {}
        OPT['alg']['type'] = 'block'
    #endif

    # Unspecified parts of M -> defaults
    if OPT['n'] >= Ny:
        raise Exception("Cannot have OPT.n larger than height of Z!")
    #endif

    if unit._utils.isempty(M):
        raise Exception("Need to specify initial model structure M!")
    else:
        M = unit.startM(Z,M)
    #endif

    # Include delays specified in model structure on inputs
    for r in range(0,nu):
        u[:,r:r+1] = np.vstack([ np.zeros([M['delay'][r], 1]) , u[0:Ny-M['delay'][r], r:r+1] ])
    #endfor

    # TODO: Pass input through a non-linearity if required by model structure
    # x = unit.u2x(u,M)
    x = u

    # Form regressor matrix
    PHI = np.empty([Ny, np.sum(M['nB']+1)])
    idx = 0
    for idu in range(0,nu):
        PHI[:, idx:idx+M['nB'][idu,0]+1] = scipy.linalg.toeplitz(x[:,idu], np.hstack([x[0:1,idu:idu+1], np.zeros([1, M['nB'][idu,0]])]))
        idx += M['nB'][idu,0]+1
    #endfor

    # Save initial model into G
    G = copy.deepcopy(M)

    # Now get the estimate via least squares (block) or some recursive method
    G['th'] = np.linalg.lstsq(PHI[OPT['n']:Ny, :], y[OPT['n']:Ny])[0]

    # Put theta into model structure
    mxB = np.max(M['nB'])
    G['B'] = np.empty([nu, mxB+1])
    idx = 0
    for idu in range(0,nu):
        G['B'][idu,:] = np.hstack([G['th'][idx:idx+M['nB'][idu,0]+1].T, np.zeros([1, mxB-M['nB'][idu,0]])])
        idx += M['nB'][idu,0]+1
    #endfor


#  % Do luxurious extras (if not doing fast version)
# if ~OPT.fast
#  % Parameter space variance of estimates:
#  pe    = y(OPT.n+1:length(y)) - PHI(OPT.n+1:length(y),:)*G.th;
#  G.var = pe'*pe/length(pe);
#  G.P   = G.var*pinv(PHI(OPT.n+1:length(y),:)'*PHI(OPT.n+1:length(y),:));
 
#  % Now load up matrix specifying standard deviations
#  P = real(sqrt(diag(G.P))); 
#  P = P(:); 
#  d = 0;
#  for r=1:nu 
#   G.SD.th(:,r) = [P(d+1:d+M.nB(r)+1);  zeros(mxB-M.nB(r),1)];
#   d = d + M.nB(r)+1;
#  end
# end


    # Load up output with model properties
    G['phi'] = PHI

    # Record that validate should use VN as the cost function to obtain prediction errors
    G['costfcn'] = 'VN'
    G['OPT'] = copy.deepcopy(OPT)

    # Add legend for prospective plotting
    G['disp'] = {}
    G['disp']['legend'] = 'Estimated n_b=' + str(mxB) + 'th order FIR model'
    G['alg'] = {}
    G['alg']['type'] = 'block' # Record that block solution was used # TODO: this should reflect actual solve type used



    return G
#endfunction

