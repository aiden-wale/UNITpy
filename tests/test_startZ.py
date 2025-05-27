
import pytest
import uonidtoolbox as unit
import numpy as np


def test_repeatedCall():
    Z = unit._utils.getExampleZData()

    Z_0 = unit.startZ(Z)
    Z_1 = pytest.matlabEng.startZ(unit._utils.data_py2ml(Z_0))

    np.testing.assert_equal(Z_0.keys(), Z_1.keys())
    np.testing.assert_equal(Z_0, Z_1)
#endfunction

@pytest.mark.matlab
def test_matlabResult():
    Z = unit._utils.getExampleZData()

    Z_py = unit.startZ(Z)
    Z_ml = pytest.matlabEng.startZ(unit._utils.data_py2ml(Z))

    np.testing.assert_equal(Z_py.keys(), Z_ml.keys())
    np.testing.assert_equal(Z_py, Z_ml)
#endfunction

