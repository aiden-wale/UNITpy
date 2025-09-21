
import uonidtoolbox as unit
import numpy as np
import scipy


def demo_fir(disp=1):

    OPT = unit.struct()
    OPT.dsp = disp

    # ====================================================
    # Specify Experiment Conditions
    # ====================================================
    T   = 0.5   # Sampling period in seconds
    var = 1e-4  # Measurement Noise Variance


    # ====================================================
    # Specify true system - ZOH impulse response of continuous defined system
    # ====================================================
    den     = np.real(np.poly([-1, -3+1j, -3-1j]))
    sysc    = scipy.signal.TransferFunction(1, den)
    sysd    = scipy.signal.cont2discrete((sysc.num, sysc.den), T, method='zoh')
    bq      = scipy.signal.dimpulse(sysd)[1][0].ravel()


    # ====================================================
    # Simulate a data record
    # ====================================================
    Z       = unit.struct()
    N       = 50*bq.shape[0]
    Z.u     = np.random.randn(1, N)
    noise   = np.sqrt(var)*np.random.randn(Z.u.size).reshape(Z.u.shape)
    Z.y     = scipy.signal.lfilter(bq, 1, Z.u) + noise


    # ====================================================
    # Specify Model Structures
    # ====================================================
    Mq          = unit.struct()
    Mq.nB       = bq.shape[0]-1
    Mq.T        = T
    Mq.type     = 'fir'


    # ====================================================
    # Estimate on basis of noise corrupted data
    # ====================================================
    Gq = unit.est(Z,Mq,OPT)


    # ====================================================
    # Plot the results
    # ====================================================
    if OPT.dsp:
        Gt              = unit.struct()
        Gt.disp         = unit.struct()
        Gt.B            = bq
        Gt.T            = T
        Gt.w            = Gq.w
        Gt.type         = 'fir'
        # Gt.colour       = 'b'
        Gt.disp.legend  = 'True Response'
        Gt.disp.legend  = Gt.disp.legend + ' q operator'
        
        unit.plotting.showbode([Gt, Gq])
    #endif

#endfunction

