
import pytest
import uonidtoolbox as unit
import numpy as np


@pytest.mark.matlab
def test_alg_fir():

    data = unit._utils.loadmat("tests/_testdata/fir.mat")
    Z   = data['Z']
    Z['u'] = Z['u'].reshape(Z['u'].size, 1)
    Z['y'] = Z['y'].reshape(Z['y'].size, 1)
    M   = data['M']
    M['in'] = np.array([M['in']])
    M['in'] = M['in'].reshape(M['in'].size)
    OPT = data['OPT']

    pytest.matlabEng.eval('load("../../tests/alg_test_data/fir_test_input.mat");', nargout=0)

    pytest.matlabEng.eval('G = fir(Z,M,OPT);', nargout=0)

    G_py = unit.fir(Z,M,OPT)
    G_ml = pytest.matlabEng.workspace['G']

    np.testing.assert_allclose(G_py['th'], G_ml['th'])

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


