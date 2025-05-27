
import numpy as np
import pytest
# import uonidtoolbox as unit


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
        dout = din.copy()
        for k in dout.keys():
            dout[k] = data_py2ml(dout[k])
        #endfor
    elif isinstance(din, np.ndarray):
        dout = din.copy()
        if 'int' in str(dout.dtype):
            dout = np.array(dout, dtype='float64')
        #endif
    elif isinstance(din, int):
        dout = float(din)
    else:
        # raise Exception("unexpected data type" + str(type(din)))
        dout = din
    #endif
    return dout
#endfunction

def helper_callMatlab_startZ(Z):
    pytest.matlabEng.workspace['Z'] = data_py2ml(Z)
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

def helper_callMatlab_startZ(Z_py):
    pytest.matlabEng.workspace['Z'] = data_py2ml(Z_py)
    pytest.matlabEng.eval('Z = startZ(Z);', nargout=0)
    Z_ml = pytest.matlabEng.workspace['Z']
    return Z_ml
#endfunction

def isempty(np_arr):
    if not isinstance(np_arr, np.ndarray):
        return False
    #endif
    if np_arr.size == 0:
        return True
    else:
        return False
    #endif
#endfunction

        