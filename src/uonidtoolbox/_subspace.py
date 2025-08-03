
import numpy as np
import uonidtoolbox as unit


def subspace(Z,M={},OPT={}):

    Z = unit.startZ(Z)

    match Z['type']:
        case 'time':
            G = unit.sid(Z,M,OPT)

        case 'frequency':
            G = unit.fsid(Z,M,OPT)

        case _:
            unit._utils.uwarning('Data type (Z.type) not known')
    #endmatch

    return G
