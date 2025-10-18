
import pytest
import uonidtoolbox as unit
import numpy as np


def test_repeatedCall():
    Z = {}
    Z['u'] = np.array([-0.3731,0.8155,0.7989,0.1202,0.5712,0.4128,-0.9870,0.7596,-0.6572,-0.6039]).reshape(1,10)
    Z['y'] = np.array([0.1769,-0.3075,-0.1318,0.5954,1.0468,-0.1980,0.3277,-0.2383,0.2296,0.4400]).reshape(1,10)

    Z_0 = unit._setup.startZ(Z)
    Z_1 = unit._setup.startZ(Z_0)

    np.testing.assert_equal(Z_0.keys(), Z_1.keys())
    np.testing.assert_equal(Z_0, Z_1)
#endfunction

@pytest.mark.matlab
def test_matlabResult():
    Z = unit.testing._utils.getExampleZData()

    Z_py = unit._setup.startZ(Z)
    Z_ml = pytest.matlabEng.startZ(unit.testing._utils.data_py2ml(Z))

    Z_py = Z_py.asdict()
    np.testing.assert_equal(Z_py.keys(), Z_ml.keys())
    np.testing.assert_equal(Z_py, Z_ml)
#endfunction

