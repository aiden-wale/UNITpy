
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
    restf   = pctrl.c2d(pctrl.tf(num,den), T, 'zoh') # introduces delay
    bq,aq   = restf.num, restf.den


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
    Gq = est(Z,Mq,OPT)


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

        mpl.figure()
        mpl.semilogx(Gt.w, )
        mpl.show()
        # showbode(Gt,Gq)
    #endif

#endfunction


