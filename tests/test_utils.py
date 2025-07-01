
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

