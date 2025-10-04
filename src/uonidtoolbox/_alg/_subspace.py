
import uonidtoolbox as unit
import numpy as np


def subspace(Z, M=unit.struct(), OPT=unit.struct()):

    Z = unit.startZ(Z)

    match Z.type:
        case 'time':
            G = unit.sid(Z,M,OPT)

        case 'frequency':
            G = unit.fsid(Z,M,OPT)

        case _:
            unit._utils.uwarning('Data type (Z.type) not known')
    #endmatch

    return G
