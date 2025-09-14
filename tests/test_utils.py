
import pytest
import uonidtoolbox as unit
import numpy as np


def test_utils_isempty_empty_dict():
    obj = {}
    res = unit._utils.isempty(obj)
    assert res == True
#endfunction

def test_utils_isempty_empty_list():
    obj = []
    res = unit._utils.isempty(obj)
    assert res == True
#endfunction

def test_utils_isempty_empty_nparray1d():
    obj = np.ndarray([0])
    res = unit._utils.isempty(obj)
    assert res == True
#endfunction

def test_utils_isempty_empty_nparray2d_0x0():
    obj = np.ndarray([0,0])
    res = unit._utils.isempty(obj)
    assert res == True
#endfunction

def test_utils_isempty_empty_nparray2d_0x1():
    obj = np.ndarray([0,1])
    res = unit._utils.isempty(obj)
    assert res == True
#endfunction

def test_utils_isempty_empty_nparray2d_1x0():
    obj = np.ndarray([1,0])
    res = unit._utils.isempty(obj)
    assert res == True
#endfunction


def test_utils_isempty_nonempty_dict():
    obj = {'a': 1}
    res = unit._utils.isempty(obj)
    assert res == False
#endfunction

def test_utils_isempty_nonempty_list():
    obj = [5]
    res = unit._utils.isempty(obj)
    assert res == False
#endfunction

def test_utils_isempty_nonempty_nparray1d():
    obj = np.ndarray([1])
    res = unit._utils.isempty(obj)
    assert res == False
#endfunction

def test_utils_isempty_nonempty_nparray2d():
    obj = np.ndarray([1,1])
    res = unit._utils.isempty(obj)
    assert res == False
#endfunction

def test_utils_isempty_nonempty_npint64():
    obj = np.int64(1)
    res = unit._utils.isempty(obj)
    assert res == False
#endfunction

def test_utils_isempty_nonempty_npfloat64():
    obj = np.float64(1)
    res = unit._utils.isempty(obj)
    assert res == False
#endfunction

def test_utils_blockhankel_2x8_4():
    u = np.empty([2,8])
    u[0,:] = np.array([1,2,3,4,5,6,7,8])
    u[1,:] = np.array([1,2,3,4,5,6,7,8])+10

    nu = u.shape[0]
    N = u.shape[1]

    nr = 4

    exp = np.empty([nu*nr, N-nr+1])
    for i in range(0, nr):
        exp[nu*i:nu*(i+1), :] = u[:, i:i+N-nr+1]
    #endfor

    res = unit._utils.blockhankel(u, nr)
    np.testing.assert_equal(res, exp)
#endfunction




