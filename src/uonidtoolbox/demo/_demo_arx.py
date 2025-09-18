
import uonidtoolbox as unit
import numpy as np
import scipy
import control as pctrl
import matplotlib.pyplot as mpl

def demo_arx(disp=1):

    OPT = unit.struct()
    OPT.dsp = disp

    # ====================================================
    # Specify Experiment Conditions
    # ====================================================
    T   = 1     # Sampling period in seconds
    N   = 100   # Number of Samples
    var = 1e-1  # Measurement Noise Variance


    # ====================================================
    # Specify true system                    
    # ====================================================
    den     = np.real(np.poly([-0.1,-1,-0.2]))
    num     = den[-1]
    restf   = pctrl.c2d(pctrl.tf(1,den), T, method='zoh') # introduces delay
    bq,aq   = restf.num[0][0], restf.den[0][0]
    bq      = np.hstack([0,bq])


    # ====================================================
    # Simulate a data record                    
    # ====================================================
    Z       = unit.struct()
    Z.u     = np.random.randn(1, N)
    noise   = np.sqrt(var)*np.random.randn(Z.u.size).reshape(Z.u.shape)
    noise   = scipy.signal.lfilter(1, aq, noise)
    Z.y     = scipy.signal.lfilter(bq, aq, Z.u) + noise


    # ====================================================
    # Specify Model Structures
    # ====================================================
    Mq          = unit.struct()
    Mq.A        = aq.shape[0]-1
    Mq.B        = bq.shape[0]-2 # delay is catered for elsewhere
    Mq.T        = T
    Mq.type     = 'arx'
    Mq.delay    = 1


    # ====================================================
    # Estimate on basis of noise corrupted data
    # ====================================================
    Gq = unit.est(Z,Mq,OPT)


    # ====================================================
    # Plot the results
    # ====================================================
    if OPT.dsp:
        Gt = unit.struct()
        Gt.A            = aq
        Gt.B            = bq
        Gt.w            = Gq.w
        Gt.type         = 'arx'
        Gt.colour       = 'b'
        Gt.disp         = unit.struct()
        Gt.disp.legend  = 'True Response'

        
        # TODO: implement all below as unit.plotting.showbode(G1, G2, ..., Gn)
        # i.e. unit.plotting.showbode(G_true, G_est)
        Gt = m2f(Gt)
        Gq = m2f(Gq)

        if Gt.G.shape != Gq.G.shape:
            raise Exception("Gt.G.shape != Gq.G.shape")

        mpl.figure()

        # Magnitude [dB]
        mpl.semilogx(Gt.w, 20*np.log10(np.abs(Gt.G)))
        mpl.semilogx(Gq.w, 20*np.log10(np.abs(Gq.G)))

        # Phase [rad/s]
        # TODO: show phase plot (as subplot)
        # TODO: provide option to show as [Hertz]

        mpl.show()
    #endif

#endfunction


# TODO: move this function out of _demo_arx
# TODO: implement handling of MIMO, MISO, SIMO systems
from numpy.polynomial.polynomial import polyval as np_polyval
import copy
def m2f(M):
    if 'finishM' not in M: M = unit.startM(M)

    G = copy.deepcopy(M)

    if M.op == 'q':
        ww      = np.exp(1j*M.w*M.T)
        pdel    = np.exp((-1j*M.w.flatten()*M.T)*G.delay) # phase lag due to delays on inputs
    else:
        raise Exception("M.op == "+str(M.op)+" not implemented")
    #endif

    match G.type:
        case "ss":
            # TODO: implement handling of SS systems
            raise Exception("m2f() not yet implemented for G.type == ss")
        case _: # polynomial model
            G.D = G.A # for M.type in ['ar','arx','arma','armax']

            # TODO: repair dimensions of A,B,C,D polynomials (unit.startM() ?)
            tmpG = unit.struct()
            for p in ['A','B','C','D']:
                if G[p].ndim > 1:
                    tmpG[p] = G[p][0,:]
                else:
                    tmpG[p] = G[p]

            A = tmpG.A[::-1]
            B = tmpG.B[::-1]

            G.G = np_polyval(1/ww, B)/np_polyval(1/ww, A)
            pp = pdel
            G.G = G.G * pp # pp is responsible for the effect of delay on input

            C = np.hstack([np.zeros(tmpG.D.size - tmpG.C.size), tmpG.C])[::-1]
            D = tmpG.D[::-1]
            G.H = np_polyval(1/ww, C)/np_polyval(1/ww, D)
    #endmatch

    # # For time series case, noise spectral factor masquerades as dynamic freq resp.
    # if nu<1: G.G = G.H

    return G
#endfunction
