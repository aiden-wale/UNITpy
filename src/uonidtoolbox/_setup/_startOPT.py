
import uonidtoolbox as unit
import numpy as np


def startOPT(OPTin=unit.struct(), Min=unit.struct()):

    o = unit.struct()

    o.n      = 0
    o.dsp    = 1
    o.miter  = 100
    o.tol    = 1e-4
    o.lmax   = 30
    o.mdec   = 1e-9
    o.fast   = 0
    o.filt   = 0
    o.step   = 1
    o.smits  = 5
    o.cost   = 'trace'
    o.subtol = np.sqrt(np.sqrt(unit._utils.eps))
    o.delta  = 1
    o.adapt  = 1
    o.dir    = 'trust'
    o.ngt    = 0
    o.nht    = 0
    o.gh     = 0
    o.saveit = 1
    o.pnum   = 100
    o.emit   = 100

    if 'type' in Min:
        if Min.type in ['ss', 'bilin', 'bilinear']:
            o.alg    = 'gn'
            o.miter  = '200'
        else:
            o.alg    = 'gn'
        #endif
    else:
        o.alg = 'gn'
    #endif
    
    o.basis  = 'polyb'
    o.smeth  = 'bilin'
    o.allP   = 1
    o.cmpgrd = 0
    o.sysnd  = 'forward'
    o.gradnd = 'midpoint'
    o.par    = 'ddlc'

    if unit._utils.isempty(OPTin):
        OPT = o.copy()
    else:
        OPT = OPTin
        for k in o.keys():
            if k not in OPT:
                OPT[k] = o[k]
            elif unit._utils.isempty(OPT[k]):
                OPT[k] = o[k]
            #endif
        #endfor
    #endif

    OPT.n = int(OPT.n)

    return OPT