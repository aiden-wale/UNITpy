
import pytest
import uonidtoolbox as unit
import numpy as np


def test_repeatedCall():
    OPT_0 = unit.startOPT()
    OPT_1 = unit.startOPT(OPT_0)

    np.testing.assert_equal(OPT_0.keys(), OPT_1.keys())
    np.testing.assert_equal(OPT_0, OPT_1)
#endfunction

@pytest.mark.matlab
def test_matlabResult():
    OPT_py = unit.startOPT()
    OPT_ml = pytest.matlabEng.startOPT()

    OPT_py = OPT_py.asdict()
    np.testing.assert_equal(OPT_py.keys(), OPT_ml.keys())
    np.testing.assert_equal(OPT_py, OPT_ml)
#endfunction
