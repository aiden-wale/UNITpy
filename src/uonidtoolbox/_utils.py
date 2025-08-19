
import numpy as np
import uonidtoolbox as unit
import warnings
import os

import sys
eps = sys.float_info.epsilon


def length(o):
    if isinstance(o, np.ndarray):
        l = np.max(o.shape)
    elif isinstance(o, int):
        l = 1
    elif isinstance(o, float):
        l = 1
    else:
        l = len(o)
    #endif
    return l
#endfunction

def isempty(obj):
    if isinstance(obj, np.ndarray):
        return True if obj.size == 0 else False
    elif isinstance(obj, (list, dict, unit.struct)):
        return True if len(obj) == 0 else False
    else:
        return False
    #endif
#endfunction


def udisp(msg, guiRunning=0, guiHandle=None):
    if guiRunning:
        #send message to GUI
        print(msg)
    else:
        print(msg)
    #endif
#endfunction

def uwarning(msg):
    warnings.warn(msg)
#endfunction

def path_to_pkg():
    # path\to\UNITpy
    return os.path.abspath(os.path.join(os.path.dirname(unit.__file__), os.pardir, os.pardir))
#endfunction


