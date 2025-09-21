
import uonidtoolbox as unit
import numpy as np
import scipy


def demo_ar(disp=1):

    OPT = unit.struct()
    OPT.dsp = disp

    # ====================================================
    # Specify Experiment Conditions
    # ====================================================
    T   = 0.5   # Sampling period in seconds
    N   = 500   # Number of Samples
    var = 1e-1  # Measurement Noise Variance


    # ====================================================
    # Specify Noise Colouring
    # ====================================================
    aq = np.real(np.poly([
        0.5, 
        0.7, 
        0.9, 
        0.9*np.exp( 1j*np.pi/6), 
        0.9*np.exp(-1j*np.pi/6)
    ]))


    # ====================================================
    # Simulate a data record
    # ====================================================
    Z       = unit.struct()
    noise   = np.sqrt(var)*np.random.randn(1,N)
    Z.y     = scipy.signal.lfilter(1, aq, noise)


    # ====================================================
    # Specify Model Structures
    # ====================================================
    Mq          = unit.struct()
    Mq.nA       = aq.shape[0]-1
    Mq.T        = T
    Mq.type     = 'ar'


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
        Gt.T            = T
        Gt.w            = Gq.w
        Gt.type         = 'ar'
        # Gt.colour       = 'b'
        Gt.disp.legend  = 'True Spectral Factor'
        
        unit.plotting.showbode([Gt, Gq])
    #endif

#endfunction

