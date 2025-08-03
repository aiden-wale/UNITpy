
import numpy as np
import scipy
import pytest
import uonidtoolbox as unit


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


def getFieldsFromMatFile(path_to_data, fieldnames):
    data = loadmat(path_to_data)
    fields = ()

    if type(fieldnames)=='str':
        fieldnames = list(fieldnames)
    #endif

    for fn in fieldnames:
        if fn not in data:
            raise Exception("data field '" + fn + "' was not found in the data")
        #endif

        match fn:
            case 'Z':
                Z = data[fn]
                Z['u'] = Z['u'].reshape(Z['u'].size, 1)
                Z['y'] = Z['y'].reshape(Z['y'].size, 1)
            case 'M':
                M = data[fn]
                M['in'] = np.array([M['in']])
                M['in'] = M['in'].reshape(M['in'].size)
            case 'OPT':
                OPT = data[fn]
            case 'G':
                G = data[fn]
                G['in'] = np.array([G['in']])
                G['in'] = G['in'].reshape(G['in'].size)
            case _:
                raise Exception("Unexpected field name. Could not retrieve data")
        #endmatch
        fields += (data[fn],)
    #endfor
        
    if len(fields) == 1:
        fields = fields[0]
    #endif

    return fields
#endfunction


# =================================================================================
# https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
# https://stackoverflow.com/posts/8832212/revisions
def loadmat(filename):
    '''
    this function should be called instead of direct scipy.io.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)
#endfunction

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], scipy.io.matlab.mat_struct):
            dict[key] = _todict(dict[key])
        #endif
    #endfor
    return dict
#endfunction       

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
        #endif
    #endfor
    return dict
#endfunction
# =================================================================================

