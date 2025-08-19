
import pytest
import uonidtoolbox as unit
import numpy as np


def test_alg_fir_fir():

    path_to_testdata = 'tests/_testdata/fir_fir.mat'

    # Z,M,OPT = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT'])
    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit.fir(Z,M,OPT)

    # pytest.matlabEng.eval('load("../../' + path_to_testdata + '");', nargout=0)
    # pytest.matlabEng.eval('G = fir(Z,M,OPT);', nargout=0)
    # G_ml = pytest.matlabEng.workspace.G

    # np.testing.assert_allclose(G_py.th, G_ml.th)
    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))

    # np.testing.assert_equal(G_py.keys(), G_ml.keys())
    # for k in G_py.keys():
    #     if k == 'in' or k == 'alg':
    #         continue
    #     #endif
    #     if isinstance(G_py[k], np.ndarray):
    #         np.testing.assert_allclose(G_py[k], G_ml[k], err_msg=k)
    #     else:
    #         np.testing.assert_equal(G_py[k], G_ml[k], err_msg=k)
    #     #endif
    # #endfor
#endfunction

def test_alg_barx_arx():

    path_to_testdata = 'tests/_testdata/barx_arx.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit.barx(Z,M,OPT)

    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
    # np.testing.assert_allclose(G_py.phi, G_ml.phi)
#endfunction

def test_alg_barx_ar():

    path_to_testdata = 'tests/_testdata/barx_ar.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit.barx(Z,M,OPT)

    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
    # np.testing.assert_allclose(G_py.phi, G_ml.phi)
#endfunction


