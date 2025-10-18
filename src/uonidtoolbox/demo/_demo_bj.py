
import uonidtoolbox as unit
import numpy as np
import scipy


def demo_bj(disp=1):

    OPT = unit.struct()
    OPT.dsp = disp

    # ====================================================
    # Specify Experiment Conditions
    # ====================================================
    T   = 1     # Sampling period in seconds
    N   = 500   # Number of Samples
    var = 1e-2  # Measurement Noise Variance


    # ====================================================
    # Specify true (linear) system
    # ====================================================
    den         = np.real(np.poly([-0.1,-1,-0.001-0.2j,-0.001+0.2j,-0.5]))
    num         = 10*den[-1]
    sysc        = scipy.signal.TransferFunction(num, den).to_ss()
    A,B,C,D,_   = scipy.signal.cont2discrete((sysc.A, sysc.B, sysc.C, sysc.D), T, method='zoh')
    bq,aq       = scipy.signal.ss2tf(A, B, C, D, 0)
    bq,aq       = np.squeeze(bq), np.squeeze(aq)
    cq          = np.array([1,-0.2])
    dq          = np.array([1,-0.5])


    # ====================================================
    # Simulate a data record
    # ====================================================
    Z       = unit.struct()
    t       = np.arange(0, N, 1)
    Z.u     = np.sign(np.sin(3*np.pi*t))
    noise   = np.sqrt(var)*np.random.randn(Z.u.size).reshape(Z.u.shape)
    noise   = scipy.signal.lfilter(cq, dq, noise)
    Z.y     = scipy.signal.lfilter(bq, aq, Z.u) + noise


    # ====================================================
    # Specify Model Structures
    # ====================================================
    Mq          = unit.struct()
    Mq.A        = aq.shape[0]-1
    Mq.B        = bq.shape[0]-2
    Mq.C        = cq.shape[0]-1
    Mq.D        = dq.shape[0]-1
    Mq.T        = T
    Mq.delay    = 1
    Mq.type     = 'bj'


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
        Gt.C            = cq
        Gt.D            = dq
        Gt.T            = T
        Gt.w            = Gq.w
        Gt.type         = 'bj'
        # Gt.colour       = 'b'
        Gt.disp.legend  = 'True Response'
        Gt.disp.legend  = Gt.disp.legend + ' q operator'
        
        unit.plotting.showbode([Gt, Gq])
    #endif

#endfunction

