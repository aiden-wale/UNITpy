
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


def test_alg_barx_ar():
    path_to_testdata = 'tests/_testdata/barx_ar.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit._alg.barx(Z,M,OPT)

    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
#endfunction


def test_alg_barx_arx():
    path_to_testdata = 'tests/_testdata/barx_arx.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit._alg.barx(Z,M,OPT)

    np.testing.assert_allclose(G_py.th, G_ml.th.reshape(G_ml.th.size,1))
#endfunction


def test_alg_gn_oe():
    path_to_testdata = 'tests/_testdata/gn_oe.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    for k in ['nA', 'nB', 'nC', 'nD']:
        if k in M: 
            if not isinstance(M[k], np.ndarray):
                M[k] = np.array(M[k])
            M[k] = M[k].reshape(M[k].size)

    OPT.dsp = 0
    G_py = unit._alg.gn(Z,M,OPT)

    np.testing.assert_equal(G_py.th.size, G_ml.th.size)
    np.testing.assert_equal(G_py.th.shape, G_ml.th.shape)

    # for k in ['A', 'B']:
    #     np.testing.assert_allclose(G_py[k], G_ml[k], rtol=1e-2)
    # #endif
    np.testing.assert_allclose(G_py.th, G_ml.th, rtol=1e-2)
#endfunction


def test_alg_gn_impl__expcurve():
    t           = np.arange(1,5,0.1)
    theta_true  = np.array([0.7394,-0.1941, 1.7119])
    y           = unit.testing._functions.expcurve(theta_true, t)
    Ny          = y.size

    def objfunc(th, compute_gradient=False):
        if not compute_gradient:
            yh = unit.testing._functions.expcurve(th, t, compute_gradient=False)
            pe = yh - y
            return 0.5*pe.ravel().dot(pe.ravel())/Ny # cost
        else:
            yh,J    = unit.testing._functions.expcurve(th, t, compute_gradient=True)
            pe      = yh - y
            cost    = 0.5*pe.ravel().dot(pe.ravel())/Ny
            g       = J.transpose() @ pe/Ny
            return cost,pe,g,J
        #endif
    #enddef

    theta0 = np.array([0.8,-0.3, 1.6])
    thetastar = unit._alg._gn._gn_impl(objfunc, theta0, maxiter=100, disp=1)

    np.testing.assert_allclose(thetastar, theta_true, rtol=1e-6)
#endfunction


def test_alg_sid_n4sid():
    path_to_testdata = 'tests/_testdata/sid_n4sid.mat'

    Z,M,OPT,G_ml = unit.testing._utils.getFieldsFromMatFile(path_to_testdata, ['Z', 'M', 'OPT', 'G'])
    G_py = unit._alg.sid(Z,M,OPT)

    G_ml.ss = unit._startM._make2d_SS_matrices(G_ml.ss, G_ml.nx, G_ml.nu, G_ml.ny)

    # Check model variables {nx, nu, ny} are equal between models
    for k in ['nx', 'nu', 'ny']: unit.testing._utils.assert_field_equal(G_py[k], G_ml[k], k, "G_py.", "G_ml.")

    # Check dimensions of system matrices
    for k in ['A', 'B', 'C', 'D']: unit.testing._utils.assert_field_equal(G_py.ss[k].shape, G_ml.ss[k].shape, k, f"G_py.ss.{k}.shape", f"G_ml.ss.{k}.shape")

    # Get transform between SS systems to disambiguate between representations
    P = unit._utils.getSimilarityTransform(G_py.ss, G_ml.ss)

    G_py.ss.A = P @ G_py.ss.A @ np.linalg.inv(P)
    G_py.ss.B = P @ G_py.ss.B
    G_py.ss.C =     G_py.ss.C @ np.linalg.inv(P)
    # G_py.ss.D = G_py.ss.D

    for k in ['A', 'B', 'C', 'D']: np.testing.assert_allclose(G_py.ss[k], G_ml.ss[k], err_msg=f"G.ss.{k}")
#endfunction

