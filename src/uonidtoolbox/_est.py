#  EST: Computes an estimate of either a Box-Jenkins model or a state space
#  model.  In the Box--Jenkins case the structure is of the form
#
#          B(p)                   C(p)
#  y_t  =  ---- u_{t-delay}   +   ---- e_t
#          A(q)                   D(p)
#
#  where e_t is white noise, p can be the forward shift operator q or
#  the Euler differencing (delta) operator d = (q-1)/T (with T being
#  the sampling period) and a quadratic (least squares) prediction error
#  loss criterion is used.  A(p) through D(p) are all polynomials in p of
#  the form
#
#  A(p) = 1.0 + a_1 p^{-1} + a_2 p^{-2} + ... + a_n p^{-n}
#  B(p) = b_0 + b_1 p^{-1} + b_2 p^{-2} + ... + b_m p^{-m}
#  C(p) = 1.0 + c_1 p^{-1} + c_2 p^{-2} + ... + c_r p^{-r}
#  D(p) = 1.0 + d_1 p^{-1} + d_2 p^{-2} + ... + d_r p^{-r}
#
#  This is the SISO linear form, but MISO forms and inclusion of
#  Hammerstein and Wiener non-linearity blocks are also supported.
#
#  In the state space case, the model is of the form
#
#  px_t = Ax_t + Bu_t + F(x_t \otimes u_t) + w_t
#   y_t = Cx_t + Du_t + G(x_t \otimes u_t) + e_t
#
#  with w_t and e_t being state and measurement noise.  In the purely
#  linear modelling case, F and G are set to zero, but otherwise a
#  bilinear model is estimated.
#
#  Usage is
#
#  G = est(Z,M,OPT);
#
#  where:
#
#   Z:          Input-Output data in one of two forms.  The standard form
#               is for it to be a record with elements Z.y and Z.u, each
#               of which are matrices with number of rows equal to the
#               number of data samples, and number of columns equal (respectively)
#               to the number of outputs and the number of inputs.  On
#               the other hand, Z can be a matrix of the form Z = [y,u]
#               where it is assumed that y is a column vector of output
#               measurements and u is a matrix whose columns are the
#               input measurements; in this latter MISO models are
#               being considered.
#
#               In Addition, both non-equidistant time domain data, and
#               frequency domain data can also be supplied via Z. Please
#               type "help startZ" at the command prompt for further
#               details.
#
#
#  M:           Data structure which defines the model structure which
#               is to be estimated from the data as follows:
#    M.A,M.B:   Initial guess for input-output dynamics.  If these are
#               given as integers, then they are interpreted as
#               specifications of the numbers of zeros in M.A and M.B.
#               Each row of M.A or M.B specifies initiation of on i/o
#               component of possible MISO model.  Any poly spec is
#               interpreted in *decreasing* powers of p, starting at
#               zeroth power.  Pad with trailing zeros if necessary to
#               make all rows of same length.
#    M.C,M.D:   Initial guess for measurement noise model.  If not
#               specified, or specified as empty matrices, then
#               the default is C/D 1; If M.C and M.D are specified as
#		             integers, then this is interpreted as the number of poles
#		             and zeros to be estimated in C and D.
#    M.delay:   Number of samples of delay to include (see above model).
#               In the case of a MISO system, this should be a vector of
#               delays, one for each input being considered.  Default is
#               delays are all equal to zero.
#    M.op:      Set to 'q' for shift and 'd' for delta.  Default = 'q'.
#    M.T:       Sampling period in s. (Ignored for q case) Default = 1;
#    M.w:       Vector of frequencies at which to calculate frequency
#               response of estimated model.  Specify in real frequency,
#               not normalised.  Default is 3 decades up to folding freq.
#    M.type:    Type of model structure for linear dynamics estimation.
#               Valid types are: 'arma','fir','arx','armax','oe','bj' for
#               polynomial models, 'ss' for a linear state space model,
#               and 'bilinear' for a bilinear state space model.
#               If M.type is not specified, it is guessed from how
#               M.A,M.B, M.C and M.D are specified.
#    M.in(k):   Structure that defines any static non-linearities on the
#               k'th input. This structure contains an element M.in(k).type
#               whose default value is `linear', but may be set to the
#               following other values in which case other elements of
#               M.in(k) (as indicated) may also be set:
#M.in(k).type:  'poly'.  This fits a polynomial model to the input
#               non-linearity, in which case M.in.eta specifies either
#               (if a single integer) the order of the polynomial or (if
#               a vector) an initial guess at the terms in the
#               polynomial.  If M.in.eta is not specified, then a default
#               cubic model is used.
#M.in(k).type:  'hinge'.  This fits a `hinging hyperplane' (ie. peicewise
#               linear) model to the input non-linearity. In which case
#               M.in.eta specifies the piecewise linearity in a fashion
#               impossible to document quickly here.  The default is a
#               deadzone shape non-linearity.
#       M.par:  Specification of state-space parametrisation used in
#               Gauss--Newton algorithm. This may be one of:
#               'ddlc' : Data-Driven-Local-Coordinates - default for
#                        state-space models (both linear and bilinear).
#               'full' : Full parametrisation.
#  OPT:         Data structure which defines options for the estimation
#               algorithm as follows:
#    OPT.alg:   Specification of estimation algorithm to use.  This may
#               be one of:
#               'gn' : Gradient based search via damped Gauss--Newton
#                      algorithm - default for SISO or MISO models.
#               'em' : Maximum likelihood estimation via EM algorithm
#                      iterative search - default for MIMO models.
#    'n4sid' or 'sid': N4SID Subspace based system identification method.
#               'cca': Canonical Correlation Analysis method.
#
#                The default is OPT.alg='gn' for SISO and MISO systems,
#                but OPT.alg='em' as a default for MIMO systems.
#    OPT.dir:   Search direction for gradient based search method, can be
#               'rgn' for robust Gauss-Newton, 'lm' for
#               Levenberg-Marquardt, 'trust' for trust region, 'bfgs' for
#               quasi-Newton method, 'grad' for pure gradient search.
#    OPT.cost:  Can be 'trace' to minimise trace of sum of squared errors,
#               or 'det' to minimise the log determinant of the sum of
#               outer product of error vectos (related to Maximum-Likelihood)
#    OPT.dsp:   Optional, set to something non-zero for verbose output.
#               Default is OPT.dsp=0;
#    OPT.n:     Number of starting data points to discard to get
#               rid of initial condition effects.  Default is none.
#    OPT.delta: Regularisation value for Levenberg-Marquardt and Trust
#               region search methods.
#    OPT.step:  Number of prediction samples ahead to use in cost
#               criterion - only applies for OPT.alg='gn'.  Default is OPT.step=1;
#    OPT.miter: Maximum number of updates of estimate from initial guess.
#               Default is 200.
#    OPT.tol:   If normalised gradient is less than OPT.tol, then gn
#               algorithm will terminate
#    OPT.lmax:  Maximum number of times search distance will be shortened
#               by bisection.  Default is 52.
#    OPT.mdec:  Minimum relative decrease of cost before search is
#               terminated.  Default is 1e-9;
#    OPT.fast:  When set to 1, then makes algorithm run maximally fast by
#               avoiding the calculation of error bounds.  Default is 0.
#
#  G:           Data structure which specifies the estimated model as
#               follows:
# G.A, G.B:     Matrices definining the estimated transfer function model.
# G.jC, G.D:     For SISO systems, these element are row vectors defining
#               co-efficients of increasing powers of M.op^-1.  For MISO,
#               they are matrices of rows, the k't row pertaining to the
#               k'th input.  For MIMO, they are 3 dim matrices with the
#               row (k,:,m) defining the transfer function from the k'th
#               input to the m'th output.
# G.ss.A,B,C:   [A,B,C,D,F,G] matrices/vectors defining estimated model in
#      D,F,G:   state space form.
# G.ss.X0,P0:   EM algorithm case only - estimate of ic's;
#    G.G:       Matrix of frequency responses.  If the system has multiple
#               inputs and multiple outputs, then this matrix is 3
#               dimensional, with one `page' per output, and the i'th
#               column of the j'th page (ie G.G(:,i,j)) being the
#               frequency response from input i to ouput j.
#    G.H:       Frequency response of estimated spectral factor of
#               measurement noise - may be multivariable as above.
#    G.Ge:      Matrix specifying 95% confidence regions for estimated
#               frequency response Ghat.  They may be plotted by using either
#               of the commands `shownyq(G)' or `showbode(G)'.
#    G.Gvar:    Vector which is var(G(w)), one element per element in G.w.
#    G.P:       Covariance Matrix for Estimated Parameters.
#    G.th:      Parameter estimates as a column vector.
#    G.mse:     Evolution of mean square cost decrease as any iterative
#               algorithm (em, gn) progresses;
#    G.LL:      EM algorithm only - evolution of log likelihood increase.
#
#   written by Brett Ninness, School of EE & CS
#              Adrian Wills   University of Newcastle
#             		          Australia.
# 
# 
# Copyright (C) Brett Ninness

import numpy as np
import uonidtoolbox as unit


def est(Z, M, OPT):

    G = 0

    # # ignore GUI for time being
    # if not Z:
    #     raise Exception("Need to specify data (Z)!")
    # elif not M:
    #     Z = unit.startZ(Z)
    #     m['nx'] = min(20, np.ceil(Z['Ny']/10))
    #     gsid = unit.subspace(Z, m)
    #     lsin = length(gsid.sing)
    #     vv = linspace(0, (lsin-1)/lsin, lsin); gsid.sing(:).T/.gsid.sing(1)
    #     [mv,mi] = min(sum(vv.*vv))
    #     M['A'] = mi + 1
    #     OPT = []
    # elif not OPT:
    #     OPT = []
    # #endif

    # # Detect if GUI running
    # gui = 0
    # guih = []
    # if 'gui' in OPT:
    #     if not isempty(OPT['gui']):
    #         gui = 1
    #         guih = OPT['gui']
    #     #endif
    # #endif

    Z   = unit.startZ(Z)
    M   = unit.startM(Z,M)
    OPT = unit.startOPT(OPT,M)
    ep  = unit.estmap(Z,M,OPT)

    if OPT['dsp']:
        dblines = "===================================================================="
        unit._utils.udisp("\n" + dblines)
        unit._utils.udisp("START ESTIMATION PROCESS:")
        unit._utils.udisp("Estimating parameters for " + M['type'].upper() + " model structure using " + M['op'].upper() + " operator.")
        unit._utils.udisp(ep['modelEquations'])
        unit._utils.udisp("INITIALISAION:")
    #endif

    # Init nonlinear parts of model if necessary
    # if not isempty(ep.startNL):
    #     M = 

    # Init estimate of system dynamics if necessary

    # Init noise model if necessary

    
    # Now call appropriate estimation algorithm

    # Fill in components of returned (estimated) model structure


    if OPT['dsp']:
        unit._utils.udisp("END ESTIMATION PROCESS")
        unit._utils.udisp(dblines + "\n")
    #endif

    return G













