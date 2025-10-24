
import uonidtoolbox as unit
import numpy as np
import scipy


def demo_arx(disp=1):

    OPT = unit._struct()
    OPT.dsp = disp

    # ====================================================
    # Specify Experiment Conditions
    # ====================================================
    T   = 1.0  # Sampling period in seconds
    N   = 300   # Number of Samples
    var = 1e-1  # Measurement Noise Variance


    # ====================================================
    # Specify true system
    # ====================================================
    den     = np.real(np.poly([-0.1,-1,-0.2]))
    num     = 1
    sysd    = scipy.signal.cont2discrete((num, den), T, method='zoh')
    bq,aq   = sysd[0].ravel(), sysd[1].ravel()


    # ====================================================
    # Simulate a data record
    # ====================================================
    Z       = unit._struct()
    Z.u     = np.random.randn(1, N)
    noise   = np.sqrt(var)*np.random.randn(Z.u.size).reshape(Z.u.shape)
    noise   = scipy.signal.lfilter(1, aq, noise)
    Z.y     = scipy.signal.lfilter(bq, aq, Z.u) + noise


    # ====================================================
    # Specify Model Structures
    # ====================================================
    Mq          = unit._struct()
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
        Gt              = unit._struct()
        Gt.disp         = unit._struct()
        Gt.A            = aq
        Gt.B            = bq
        Gt.w            = Gq.w
        Gt.type         = 'arx'
        # Gt.colour       = 'b'
        Gt.disp.legend  = 'True Response'

        unit.plotting.showbode([Gt, Gq])
    #endif

#endfunction

