
import numpy as np
import pytest
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
    else:
        l = len(o)
    #endif
    return l
#endfunction

def getExampleZData():
    Z = np.arange(4*6).reshape(4,6)
    return Z
#endfunction


def data_py2ml(din):
    if isinstance(din, dict):
        for k in din.keys():
            din[k] = data_py2ml(din[k])
        #endfor
    elif isinstance(din, np.ndarray):
        if 'int' in str(din.dtype): # accounts for any type with 'int', i.e. int32, int64...
            din = np.array(din, dtype='float64')
        #endif
    elif isinstance(din, int):
        din = float(din)
    else:
        # din = din
        pass
    #endif
    return din
#endfunction

# def data_ml2py(din):
#     if isinstance(din, dict):
#         for k in din.keys():
#             din[k] = data_ml2py(din[k])
#         #endfor
#     elif isinstance(din, matlab.double):
#         din = np.array(din, dtype='float64')
#     else:
#         # raise Exception("unexpected data type" + str(type(din)))
#         din = din
#     #endif
#     return din
# #endfunction


def helper_callMatlab_startZ(Z_py):
    pytest.matlabEng.workspace['Z'] = data_py2ml(Z_py)
    pytest.matlabEng.eval('Z = startZ(Z);', nargout=0)
    Z_ml = pytest.matlabEng.workspace['Z']
    return Z_ml
#endfunction

def helper_callMatlab_startM(*args):
    if len(args) == 0:
        pytest.matlabEng.eval('M = startM();', nargout=0)
    elif len(args) == 2:
        pytest.matlabEng.workspace['Z'] = data_py2ml(args[0])
        pytest.matlabEng.workspace['M'] = data_py2ml(args[1])
        pytest.matlabEng.eval('M = unitpy_test_helper_startM_py2ml(M);', nargout=0)
        pytest.matlabEng.eval('M = startM(Z,M);', nargout=0)
    else:
        raise Exception("helper_callMatlab_startM() must have 0 or 2 arguments(Z,M)")
    #endif
    pytest.matlabEng.eval('M = unitpy_test_helper_startM_ml2py(M);', nargout=0)
    M_ml = pytest.matlabEng.workspace['M']
    return M_ml
#endfunction

def helper_callMatlab_estmap(Z,M,OPT):
    pytest.matlabEng.workspace['Z'] = data_py2ml(Z)
    pytest.matlabEng.workspace['M'] = data_py2ml(M)
    pytest.matlabEng.workspace['OPT'] = data_py2ml(OPT)
    pytest.matlabEng.eval('M = unitpy_test_helper_startM_py2ml(M);', nargout=0)
    pytest.matlabEng.eval('ep = estmap(Z,M,OPT);', nargout=0)
    ep_ml = pytest.matlabEng.workspace['ep']
    return ep_ml
#endfunction

def isempty(obj):
    if isinstance(obj, np.ndarray):
        return True if obj.size == 0 else False
    elif isinstance(obj, (list, dict)):
        return True if len(obj) == 0 else False
    else:
        return False
    #endif
#endfunction


def udisp(msg, guiRunning=0, guiHandle=[]):
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



import scipy

def loadmat(filename):
    '''
    this function should be called instead of direct scipy.io.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], scipy.io.matlab.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, scipy.io.matlab.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

