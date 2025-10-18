
import uonidtoolbox as unit
import numpy as np
import scipy
import copy


def VN(theta, Z, M, OPT, compute_gradient=False):
    # TODO: handle MISO polynomial case(s)
    # Extract inputs and outputs specified
    y,u,ny,nu,Ny = unit._setup._startZ._Z2data(Z)

    # Include delays specified in model structure on inputs
    for r in range(0,nu):
        u[:,r:r+1] = np.vstack([ np.zeros([M.delay[r], 1]) , u[0:Ny-M.delay[r], r:r+1] ])
    #endfor

    # Get polynomials in model structure form, from theta
    Mn = unit._utils.theta2m(theta, M)
    a = Mn.A.ravel()
    b = Mn.B.ravel()
    c = Mn.C.ravel()
    d = Mn.D.ravel()

    yh = scipy.signal.lfilter(b, a, u.transpose()).transpose()

    # Compute prediction errors and cost
    pe      = yh - y
    cost    = pe.ravel().dot(pe.ravel())/Ny

    if not compute_gradient:
        return cost
    #endif

    # Compute jacobian
    J = np.ndarray([Ny, theta.size])
    idx = 0
    ac = np.convolve(a, c)
    for k in range(0, M.nB[0]+1):
        num = np.zeros(M.nB[0]+1); num[k] = 1
        J[:, idx] = scipy.signal.lfilter(np.convolve(num,d), ac, u.transpose()).ravel()
        idx += 1
    #endfor
    for k in range(1, M.nA[0]+1):
        num = np.zeros(M.nA[0]+1); num[k] = 1
        J[:, idx] = -scipy.signal.lfilter(np.convolve(d,num), ac, yh.transpose()).ravel()
        idx += 1
    #endfor
    if M.type == 'bj':
        for k in range(1, M.nC[0]+1):
            num = np.zeros(M.nC[0]+1); num[k] = 1
            J[:, idx] = scipy.signal.lfilter(num, c, pe.transpose()).ravel()
            idx += 1
        #endfor
        for k in range(1, M.nD[0]+1):
            num = np.zeros(M.nD[0]+1); num[k] = 1
            J[:, idx] = -scipy.signal.lfilter(num, ac, yh.transpose()).ravel()
            idx += 1
        #endfor
    #endif

    g = 2*J.transpose() @ pe/Ny

    return cost,pe,g,J
#enddef


class IterationInfoPrinter:
    t = []
    w = []

    def __init__(self, col_titles=[], min_col_widths=[]):
        if len(col_titles) < 1:
            raise Exception("Must have titles for each column.")

        if len(min_col_widths) > 0: # If minimum col widths have been provided...
            if len(min_col_widths) != len(col_titles): # but doesn't match number of titles
                raise Exception("Number of titles must be same number of provided column widths.")

        self._setColumnTitles(col_titles)
        self._setColumnWidths(min_col_widths)
    #enddef

    def _setColumnTitles(self, col_titles):
        self.t = col_titles
    #enddef

    def _setColumnWidths(self, min_col_widths):
        nc = len(self.t)
        self.w = [0]*nc if (len(min_col_widths) != nc) else min_col_widths

        for i in range(0, nc):
            if self.w[i] < len(self.t[i]) + 2:
                self.w[i] = len(self.t[i]) + 2
            #endif
        #endfor
    #enddef

    def printHeader(self):
        self._printLine('top')

        headerstr = f" {self.t[0]:>{self.w[0]}}"
        for i in range(1, len(self.w)):
            headerstr += f" \u2502 {self.t[i]:>{self.w[i]}}"
        #endfor
        unit._utils.udisp(headerstr)

        self._printLine('mid')
    #enddef

    def _printLine(self, linepos):
        lsepr = "\u252c" if linepos=='top' else "\u256a" if linepos=='mid' else "\u2534" # if 'btm'
        lflat = "\u2500" if linepos=='top' else "\u2550" if linepos=='mid' else "\u2500" # if 'btm'
        linestr = lflat*(self.w[0]+2)
        for i in range(1, len(self.w)):
            linestr += lsepr + lflat*(self.w[i]+2)
        #endfor
        unit._utils.udisp(linestr)
    #enddef

    def printDataRow(self, data):
        d = data
        infostr = f" {d[0]:>{self.w[0]}}"
        for i in range(1, len(d)):
            infostr += f" \u2502 {d[i]:>{self.w[i]}.5e}" if isinstance(d[i], (int, float)) else f" \u2502 {d[i]:>{self.w[i]}}"
        #endfor
        unit._utils.udisp(infostr)
    #enddef

    def finishTable(self):
        self._printLine('btm')
    #enddef
#endclass
