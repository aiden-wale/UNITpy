
import numpy as np
import uonidtoolbox as unit

length = unit._utils.length
isempty = unit._utils.isempty


def startNL(Z, M):

    if 'in' not in M:
        M['in'] = []
        for k in range(0, M['nu']):
            M['in'].append({})
            M['in'][k]['type'] = 'linear'
        #endfor
    elif isempty(M['in']):
        for k in range(0, M['nu']):
            M['in'].append({})
            M['in'][k]['type'] = 'linear'
        #endfor
    else:
        for k in range(0, M['nu']):
            if 'type' not in M['in'][k]:
                M['in'][k]['type'] = 'linear'
            #endif
        #endfor
    #endif
    if 'out' not in M:
        M['out'] = {}
        M['out']['type'] = 'linear'
    elif 'type' not in M['out']:
        M['out']['type'] = 'linear'
    #endif

    if M['nu'] == 0:
        M['in'] = np.array([]).reshape([0])
    #endif
    if M['ny'] == 0:
        M['out'] = np.array([]).reshape([0])
    #endif

    # Loop over input nonlinearities to specify default values
    for k in range(0, M['nu']):
        match M['in'][k]['type']:
            case 'linear':
                M['in'][k]['eta']   = np.array([]).reshape([0,0])
                M['in'][k]['neta']  = 0

            case 'hinge':
                # If hinging hyperplane non-linearity, then remove spurious specs for other types.
                if 'upper' in M['in'][k]: M['in'][k]['upper'] = np.array([]).reshape([0,0])
                if 'lower' in M['in'][k]: M['in'][k]['lower'] = np.array([]).reshape([0,0])
                if ('eta' not in M['in'][k]) or ('eta' in M['in'][k] and isempty(M['in'][k]['eta'])):
                    M['in'][k]['eta'] = np.array([0.05, 1, -0.05, -1, -0.05, 1])
                else:
                    if length(M['in'][k]['eta']) % 2:
                        raise Exception("Initial guess in M.in.eta must contain pairs of parameters and hence must be of even length")
                    #endif
                #endif
                M['in'][k]['neta'] = length(M['in'][k]['eta'])

            case 'poly':
                # If polynomial non-linearity, then remove spurious specs for other types.
                if 'upper' in M['in'][k]: M['in'][k]['upper'] = np.array([]).reshape([0,0])
                if 'lower' in M['in'][k]: M['in'][k]['lower'] = np.array([]).reshape([0,0])
                if ('eta' not in M['in'][k]) or ('eta' in M['in'][k] and isempty(M['in'][k]['eta'])):
                    M['in'][k]['eta'] = np.vstack([ [[1]], np.zeros([11,1]) ])
                elif length(M['in'][k]['eta']) == 1:
                    M['in'][k]['eta'] = np.vstack([ [[1]], np.zeros([int(np.floor(M['in'][k]['eta'])-1), 1]) ])
                #endif
                M['in'][k]['neta'] = length(M['in'][k]['eta'])

            case ('saturation', 'deadzone'):
                if 'eta' in M['in'][k]:         M['in'][k]['eta']   = np.array([]).reshape([0,0])
                if 'upper' not in M['in'][k]:   M['in'][k]['upper'] = np.array([]).reshape([0,0])
                if 'lower' not in M['in'][k]:   M['in'][k]['lower'] = np.array([]).reshape([0,0])
                if isempty(M['in'][k]['lower']) and isempty(M['in'][k]['upper']):
                    M['in'][k]['upper'] =  0.5
                    M['in'][k]['lower'] = -0.5
                    M['in'][k]['neta']  =  2; 
                elif not isempty(M['in'][k]['lower']) and not isempty(M['in'][k]['upper']):
                    # Sanity check its specs for saturation and deadzone non-linearities
                    if M['in'][k]['lower'] > M['in'][k]['upper']:
                        raise Exception("Must have M.in.lower < M.in.upper on input #" + str(k))
                    #endif
                    M['in'][k]['neta'] = 2
                else:
                    if not isempty(M['in'][k]['lower']):
                        M['in'][k]['upper'] = -M['in'][k]['lower']
                        M['in'][k]['neta']  = 1
                    else:
                        M['in'][k]['lower'] = -M['in'][k]['upper']
                        M['in'][k]['neta']  = 1
                    #endif
                #endif

            case _:
                raise Exception("Error: Input non-linearity on input " + str(k) + " must be one of linear, saturation, deadzone, hinge or poly")

        #endmatch   
    #endfor

    # Now set defaults for output nonlinearity
    match M['out']['type']:
        case 'linear':
            M['out']['eta']     = np.array([]).reshape([0,0])
            M['out']['neta']    = 0

        case 'hinge':
            # If hinging hyperplane non-linearity, then remove spurious specs for other types.
            if 'upper' in M['out']: del M['out']['upper']
            if 'lower' in M['out']: del M['out']['lower']
            if 'eta' not in M['out']:
                M['out']['eta'] = np.array([0, 1, 0, 0, 0, 0])
            #endif
            M['out']['neta'] = length(M['out']['eta'])

        case 'poly':
            # If polynomial non-linearity, then remove spurious specs for other types.
            if 'upper' in M['out']: del M['out']['upper']
            if 'lower' in M['out']: del M['out']['lower']
            if 'eta' not in M['out']:
                M['out']['eta'] = np.array([0.1, 1, 0.1, 0.1, 0.01, 0.01])
            elif length(M['out']['eta']) == 1:
                # M['out']['eta'][k] = np.concat([ [[1]], np.zeros([int(np.floor(M['out']['eta'])-1), 1]) ]) # TODO: is this right? from startNL.m, uses 'k' index despite 'k' out of scope
                M['out']['eta'] = np.vstack([ [[1]], np.zeros([int(np.floor(M['out']['eta'])-1), 1]) ])
            #endif
            M['out']['neta'] = length(M['out']['eta'])

        case ('saturation', 'deadzone'):
            if 'eta' in M['out']: del M['out']['eta']
            if 'lower' not in M['out'] and 'upper' not in M['out']:
                M['out']['upper'] =  0.1
                M['out']['lower'] = -0.1
                M['out']['neta']  =  2; 
            elif 'lower' in M['out'] and 'upper' in M['out']:
                # Sanity check its specs for saturation and deadzone non-linearities
                if M['out']['lower'] > M['out']['upper']:
                    raise Exception("Must have M.out.lower < M.out.upper on input")
                #endif
                M['out']['neta'] = 2
            else:
                if 'lower' in M['out']:
                    M['out']['upper'] = -M['out']['lower']
                    M['out']['neta']  = 1
                else:
                    M['out']['lower'] = -M['out']['upper']
                    M['out']['neta']  = 1
                #endif
            #endif

        case _:
            raise Exception("Error: Output non-linearity must be one of linear, saturation, deadzone, hinge or poly")

    #endmatch  
    
    return M
    