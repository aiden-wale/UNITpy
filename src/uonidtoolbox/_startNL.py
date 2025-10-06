
import uonidtoolbox as unit
import numpy as np

length = unit._utils.length
isempty = unit._utils.isempty


def startNL(Z, M):

    if 'inp' not in M:
        M.inp = []
        for k in range(0, M.nu):
            M.inp.append(unit.struct())
            M.inp[k].type = 'linear'
        #endfor
    elif isempty(M.inp):
        for k in range(0, M.nu):
            M.inp.append(unit.struct())
            M.inp[k].type = 'linear'
        #endfor
    else:
        for k in range(0, M.nu):
            if 'type' not in M.inp[k]:
                M.inp[k].type = 'linear'
            #endif
        #endfor
    #endif
    if 'out' not in M:
        M.out = unit.struct()
        M.out.type = 'linear'
    elif 'type' not in M.out:
        M.out.type = 'linear'
    #endif

    if M.nu == 0:
        M.inp = np.array([]).reshape([0]) #TODO: should be list of size 0, i.e. M.inp = []
    #endif
    if M.ny == 0:
        M.out = np.array([]).reshape([0])
    #endif

    # Loop over input nonlinearities to specify default values
    for k in range(0, M.nu):
        match M.inp[k].type:
            case 'linear':
                M.inp[k].eta   = np.array([]).reshape([0,0])
                M.inp[k].neta  = 0

            case 'hinge':
                # If hinging hyperplane non-linearity, then remove spurious specs for other types.
                if 'upper' in M.inp[k]: M.inp[k].upper = np.array([]).reshape([0,0])
                if 'lower' in M.inp[k]: M.inp[k].lower = np.array([]).reshape([0,0])
                if ('eta' not in M.inp[k]) or ('eta' in M.inp[k] and isempty(M.inp[k].eta)):
                    M.inp[k].eta = np.array([0.05, 1, -0.05, -1, -0.05, 1])
                else:
                    if length(M.inp[k].eta) % 2:
                        raise Exception("Initial guess in M.inp.eta must contain pairs of parameters and hence must be of even length")
                    #endif
                #endif
                M.inp[k].neta = length(M.inp[k].eta)

            case 'poly':
                # If polynomial non-linearity, then remove spurious specs for other types.
                if 'upper' in M.inp[k]: M.inp[k].upper = np.array([]).reshape([0,0])
                if 'lower' in M.inp[k]: M.inp[k].lower = np.array([]).reshape([0,0])
                if ('eta' not in M.inp[k]) or ('eta' in M.inp[k] and isempty(M.inp[k].eta)):
                    M.inp[k].eta = np.vstack([ [[1]], np.zeros([11,1]) ])
                elif length(M.inp[k].eta) == 1:
                    M.inp[k].eta = np.vstack([ [[1]], np.zeros([int(np.floor(M.inp[k].eta)-1), 1]) ])
                #endif
                M.inp[k].neta = length(M.inp[k].eta)

            case ('saturation', 'deadzone'):
                if 'eta' in M.inp[k]:         M.inp[k].eta   = np.array([]).reshape([0,0])
                if 'upper' not in M.inp[k]:   M.inp[k].upper = np.array([]).reshape([0,0])
                if 'lower' not in M.inp[k]:   M.inp[k].lower = np.array([]).reshape([0,0])
                if isempty(M.inp[k].lower) and isempty(M.inp[k].upper):
                    M.inp[k].upper =  0.5
                    M.inp[k].lower = -0.5
                    M.inp[k].neta  =  2; 
                elif not isempty(M.inp[k].lower) and not isempty(M.inp[k].upper):
                    # Sanity check its specs for saturation and deadzone non-linearities
                    if M.inp[k].lower > M.inp[k].upper:
                        raise Exception("Must have M.inp.lower < M.inp.upper on input #" + str(k))
                    #endif
                    M.inp[k].neta = 2
                else:
                    if not isempty(M.inp[k].lower):
                        M.inp[k].upper = -M.inp[k].lower
                        M.inp[k].neta  = 1
                    else:
                        M.inp[k].lower = -M.inp[k].upper
                        M.inp[k].neta  = 1
                    #endif
                #endif

            case _:
                raise Exception("Error: Input non-linearity on input " + str(k) + " must be one of linear, saturation, deadzone, hinge or poly")

        #endmatch   
    #endfor

    # Now set defaults for output nonlinearity
    match M.out.type:
        case 'linear':
            M.out.eta     = np.array([]).reshape([0,0])
            M.out.neta    = 0

        case 'hinge':
            # If hinging hyperplane non-linearity, then remove spurious specs for other types.
            if 'upper' in M.out: del M.out.upper
            if 'lower' in M.out: del M.out.lower
            if 'eta' not in M.out:
                M.out.eta = np.array([0, 1, 0, 0, 0, 0])
            #endif
            M.out.neta = length(M.out.eta)

        case 'poly':
            # If polynomial non-linearity, then remove spurious specs for other types.
            if 'upper' in M.out: del M.out.upper
            if 'lower' in M.out: del M.out.lower
            if 'eta' not in M.out:
                M.out.eta = np.array([0.1, 1, 0.1, 0.1, 0.01, 0.01])
            elif length(M.out.eta) == 1:
                # M.out.eta[k] = np.concat([ [[1]], np.zeros([int(np.floor(M.out.eta)-1), 1]) ]) # TODO: is this right? from startNL.m, uses 'k' index despite 'k' out of scope
                M.out.eta = np.vstack([ [[1]], np.zeros([int(np.floor(M.out.eta)-1), 1]) ])
            #endif
            M.out.neta = length(M.out.eta)

        case ('saturation', 'deadzone'):
            if 'eta' in M.out: del M.out.eta
            if 'lower' not in M.out and 'upper' not in M.out:
                M.out.upper =  0.1
                M.out.lower = -0.1
                M.out.neta  =  2; 
            elif 'lower' in M.out and 'upper' in M.out:
                # Sanity check its specs for saturation and deadzone non-linearities
                if M.out.lower > M.out.upper:
                    raise Exception("Must have M.out.lower < M.out.upper on input")
                #endif
                M.out.neta = 2
            else:
                if 'lower' in M.out:
                    M.out.upper = -M.out.lower
                    M.out.neta  = 1
                else:
                    M.out.lower = -M.out.upper
                    M.out.neta  = 1
                #endif
            #endif

        case _:
            raise Exception("Error: Output non-linearity must be one of linear, saturation, deadzone, hinge or poly")

    #endmatch  
    
    return M
    