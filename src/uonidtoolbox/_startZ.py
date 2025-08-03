
import numpy as np
import uonidtoolbox as unit


def startZ(Z):

    for k in ['Ny','ny','nu','passed_startZ']:
        if k in Z:
            Z[k] = int(Z[k])
        #endif
    #endfor

    if 'passed_startZ' in Z:
        return Z

    if not isinstance(Z, dict) and not isinstance(Z, np.ndarray):
        raise Exception("Z must be a dictionary or numpy.ndarray")

    if isinstance(Z, dict):
        if 'y' not in Z:
            raise Exception("Must have a Z.y field")
    else:
        (N, nin) = Z.shape
        if N<nin:
            Z = Z.T
            (N, nin) = Z.shape

        Zm = Z
        Z = {}
        y = Zm[:, 0]
        y = y.reshape(y.size,1)
        if nin>1:
            if np.all(np.isreal(y)):
                Z['y'] = y
                Z['u'] = Zm[:,1:nin]
            else:
                Z['y'] = y.reshape(1,1,y.size)
                Z['w'] = Zm[:,1]
        else:
            if np.all(np.isreal(y)):
                Z['y'] = y
            else:
                raise Exception("First column of data is complex valued, but there is no second column (frequency).")

    if np.any(Z['y'].shape == 0):
        raise Exception("One or more dimensions of Z.y is 0. Must have data to proceed.")

    if 'type' not in Z:
        if np.all(np.isreal(Z['y'])):
            Z['type'] = 'time'
        else:
            Z['type'] = 'frequency'
    else:
        match Z['type']:
            case 'time':
                if not np.all(np.isreal(Z['y'])):
                    raise Exception("startZ:typeInconsistentReal', 'Z.type is inconsistent with Z.y (Z.y should be real valued).")
            case 'frequency':
                if np.all(np.isreal(Z['y'])):
                    raise Exception("startZ:typeInconsistentComplex', 'Z.type is inconsistent with Z.y (Z.y should be complex valued).")
            case _:
                raise Exception("The value in Z.type is not known.")
    
    match Z['type']:
        # --------------------------------------------------------------------------
        #    TIME DOMAIN DATA
        # --------------------------------------------------------------------------
        case 'time':
            (Ny,p) = Z['y'].shape
            if p>Ny: Z['y'] = Z['y'].T
            (Ny,p) = Z['y'].shape
            Z['Ny'] = Ny
            Z['ny'] = p
            Z['nu'] = 0

            if 'u' in Z:
                (Nu,m) = Z['u'].shape
                if np.any(Z['y'].shape == 0):
                    Z['u'] = np.array([])
                else:
                    if m>Nu: Z['u'] = Z['u'].T
                    (Nu,m) = Z['u'].shape
                    Z['nu'] = m
                    if Nu != Ny:
                        raise Exception("Number of input samples does NOT equal number of output samples")
            else:
                Z['u'] = np.array([])

            if 'T' in Z:
                if Z['T'] > 0: # data should be discrete
                    if 't' in Z:
                        # make sure t is equidistant
                        if np.abs(np.max(np.diff(Z['t'])) - np.min(np.diff(Z['t']))) > 1e-10:
                            raise Exception("Z.t is not equidistant, but Z.T>0")
                    else:
                        Z['t'] = np.arange(0, Ny*Z['T'], Z['T']).reshape(1,Ny)
                else: # data should be continuous
                    if 't' not in Z:
                        raise Exception("Time stamps Z.t must be supplied for non-equidistant data (i.e. if Z.T=0)")

            else: # no Z.T field, so try and decide
                if 't' in Z:
                    if np.max(np.diff(Z['t'])) - np.min(np.diff(Z['t'])) > 1e-12:
                        Z['T'] = 0
                    else:
                        Z['T'] = np.mean(np.diff(Z['t']))
                else:
                    Z['T'] = 1
                    Z['t'] = np.arange(0, Ny*Z['T'], Z['T']).reshape(1,Ny)

            if 'D' not in Z:
                Z['D'] = Z['T']

            if 'd' in Z:
                if np.max(Z['d']) > np.min(np.diff(Z['t'])):
                    raise Exception("An integration time in Z.d is greater than the time between samples in Z.t!")
                elif Z['d'].size > 1 and Z['d'].size < N:
                    raise Exception("Must have the same number of entries in integration times Z.d as there are samples!")
                elif Z['d'].size == 1:
                    Z['d'] = Z['d']*np.ones([1,Z['Ny']])
            else:
                Z['d'] = np.min(np.diff(Z['t']))*np.ones([1,Z['Ny']])

        # --------------------------------------------------------------------------
        #    FREQUENCY DOMAIN DATA
        # --------------------------------------------------------------------------
        case 'frequency':
            if Z['y'].ndim < 3:
                (M, N) = Z['y'].shape
                if M>N:
                    Z['y'] = Z['y'].T
                (M, N) = Z['y'].shape
                if M<2:
                    Z['y'] = Z['y'].reshape(1,1,Z['y'].size)
                else:
                    raise Exception("startZ:lackOfDimensions', 'Z.y does not have enough  dimensions")

            (p,m,Ny) = Z['y'].shape
            if p>Ny or m>Ny:
                raise Exception("Z.y needs to be transposed")
            Z['Ny'] = Ny
            Z['ny'] = p
            Z['nu'] = m

            if 'u' in Z:
                if Z['u'].ndim < 3:
                    (M, N) = Z['u'].shape
                    if M>N:
                        Z['u'] = Z['u'].T
                    (M, N) = Z['u'].shape
                    if M<2:
                        Z['u'] = Z['u'].reshape(1,1,Z['u'].size)
                    else:
                        # Assume a multi-input signal that is given as a m x N matrix
                        Z['u'] = Z['u'].reshape(M,1,N)
            else:
                Z['u'] = np.zeros(Z['ny'], Z['nu'], Z['Ny'])
                for k in range(0, Z['Ny']):
                    Z['u'][:,:,k] = np.eye(Z['ny'], Z['nu'])

            if 'w' not in Z:
                raise Exception("Frequency domain data requires a Z.w field")
            else:
                if np.max(np.abs(Z['w'])) > np.pi:
                    if 'T' in Z:
                        if np.max(np.abs(Z['w'])) > np.pi/Z['T'] + np.sqrt(unit._utils.eps):
                            raise Exception("Frequency range is greater than expected, i.e. max(abs(Z.w)) > pi/Z.T.")
                    else: 
                        Z['T'] = 0
                else: # should be discrete
                    if 'T' not in Z:
                        Z['T'] = 1

        case _:
            raise Exception("Value in Z.type is not known")

    Z['passed_startZ'] = True


    return Z;




def Z2data(Z):
    Z = startZ(Z)

    match Z['type']:
        case 'time':
            y   = Z['y']
            u   = Z['u']
            ny  = Z['ny']
            nu  = Z['nu']
            Ny  = Z['Ny']

        case 'frequency':
            y   = Z['y']
            u   = Z['w']
            ny  = Z['ny']
            nu  = Z['nu']
            Ny  = Z['Ny']

        case _:
            raise Exception("Value in Z.type not known.")
    #endmatch

    return y,u,ny,nu,Ny,Z
#endfunction