
import pytest
import uonidtoolbox as unit
import numpy as np


def test_alg_fir_fir():
    path_to_testdata = 'tests/_testdata/fir_fir.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit.fir(Z,M,OPT)

    # np.testing.assert_allclose(G_py.th, G_ml.th)
    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
#endfunction

def test_alg_barx_arx():
    path_to_testdata = 'tests/_testdata/barx_arx.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit.barx(Z,M,OPT)

    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
#endfunction

def test_alg_barx_ar():
    path_to_testdata = 'tests/_testdata/barx_ar.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit.barx(Z,M,OPT)

    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
#endfunction

@pytest.mark.skip
def test_alg_sid():
    assert False
#endfunction

@pytest.mark.skip
def test_alg_fsid():
    assert False
#endfunction