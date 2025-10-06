
import uonidtoolbox as unit
import numpy as np
import scipy


def demo_oe(disp=1):

    OPT = unit.struct()
    OPT.dsp = disp

    # ====================================================
    # Specify Experiment Conditions
    # ====================================================
    T   = 1   # Sampling period in seconds
    N   = 500   # Number of Samples
    var = 1e-1  # Measurement Noise Variance


    # ====================================================
    # Specify true (linear) system
    # ====================================================
    den         = np.real(np.poly([-0.1,-1,-0.2,-0.3,-0.5,-0.05+1j*3,-0.05-1j*3]))
    num         = 10*den[-1]
    sysc        = scipy.signal.TransferFunction(num, den).to_ss()
    A,B,C,D,_   = scipy.signal.cont2discrete((sysc.A, sysc.B, sysc.C, sysc.D), T, method='zoh')
    bq,aq       = scipy.signal.ss2tf(A, B, C, D, 0)
    bq,aq       = np.squeeze(bq), np.squeeze(aq)


    # ====================================================
    # Simulate a data record
    # ====================================================
    Z       = unit.struct()
    Z.u     = np.sign(np.sin(3*np.pi*np.arange(0, N-1, 1)))
    noise   = np.sqrt(var)*np.random.randn(Z.u.size).reshape(Z.u.shape)
    Z.y     = scipy.signal.lfilter(bq, aq, Z.u) + noise


    # ====================================================
    # Specify Model Structures
    # ====================================================
    Mq          = unit.struct()
    Mq.A        = aq.shape[0]-1
    Mq.B        = Mq.A-1
    Mq.T        = T
    Mq.delay    = 1
    Mq.type     = 'oe'


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
        Gt.A            = aq
        Gt.B            = bq
        Gt.T            = T
        Gt.w            = Gq.w
        Gt.type         = 'oe'
        # Gt.colour       = 'b'
        Gt.disp.legend  = 'True Response'
        Gt.disp.legend  = Gt.disp.legend + ' q operator'
        
        unit.plotting.showbode([Gt, Gq])
    #endif

#endfunction

