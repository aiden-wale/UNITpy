
import pytest
import uonidtoolbox as unit
import numpy as np


def test_alg_fir_fir():
    path_to_testdata = 'tests/_testdata/fir_fir.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit._alg.fir(Z,M,OPT)

    # np.testing.assert_allclose(G_py.th, G_ml.th)
    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
#endfunction

def test_alg_barx_arx():
    path_to_testdata = 'tests/_testdata/barx_arx.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit._alg.barx(Z,M,OPT)

    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
#endfunction

def test_alg_barx_ar():
    path_to_testdata = 'tests/_testdata/barx_ar.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit._alg.barx(Z,M,OPT)

    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
#endfunction

def test_alg_sid():
    path_to_testdata = 'tests/_testdata/sid_n4sid.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit._alg.sid(Z,M,OPT)

    G_ml.ss = unit._startM._make2d_SS_matrices(G_ml.ss, G_ml.nx, G_ml.nu, G_ml.ny)

    # Check model variables {nx, nu, ny} are equal between models
    for k in ['nx', 'nu', 'ny']: unit.testing._utils.assert_field_equal(G_py[k], G_ml[k], k, "G_py.", "G_ml.")

    # Check dimensions of system matrices
    for k in ['A', 'B', 'C', 'D']: unit.testing._utils.assert_field_equal(G_py.ss[k].shape, G_ml.ss[k].shape, k, "G_py.ss."+k+".shape", "G_ml.ss."+k+".shape")

    # Get transform between SS systems to disambiguate between representations
    T = unit._utils.getSimilarityTransform(G_py.ss, G_ml.ss)

    G_py.ss.A = T @ G_py.ss.A @ np.linalg.inv(T)
    G_py.ss.B = T @ G_py.ss.B
    G_py.ss.C =     G_py.ss.C @ np.linalg.inv(T)
    # G_py.ss.D = G_py.ss.D

    for k in ['A', 'B', 'C', 'D']: np.testing.assert_allclose(G_py.ss[k], G_ml.ss[k], err_msg="G.ss."+k)
#endfunction

@pytest.mark.skip
def test_alg_fsid():
    assert False
#endfunction