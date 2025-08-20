
import numpy as np
import uonidtoolbox as unit

isempty = unit._utils.isempty


def estmap(Z, M, OPT):

    ep = unit.struct()

    ep.startG    = ""
    ep.startH    = ""
    ep.startNL   = ""
    ep.alg       = ""
    ep.finishM   = "finishM"

    if isempty(OPT):
        unit._utils.uwarning("estmap: 3 inputs should be supplied")
        return ep
    #endif

    if isempty(M):
        return ep
    elif not isinstance(M, (dict, unit.struct)):
        return ep
    elif 'type' not in M:
        return ep
    #endif

    # Setup strings so that we can print model equations
    if 'inp' in M:
        nli = 0
        for i in range(0, M.nu):
            if not M.inp[i].type == 'linear':
                nli = 1
            #endif
        #endfor
        if nli:
            ham = ''
            if M.nu > 1:
                inp = 'x_i(t)'
                for i in range(0, M.nu):
                    ham = ham + '   x_%i(t) = %s(u_%i(t))\n' %(i, M.inp[i].type, i)
                #endfor
            else:
                inp = 'x(t)'
                ham = '     x(t) = %s(u(t))\n' %(M.inp[i].type)
            #endif
        else:
            ham = ''
            if M.nu > 1:
                inp = 'u_t(t)'
            else:
                inp = 'u(t)'
            #endif
        #endif
    #endif

    if 'out' in M:
        if not M.out.type == 'linear':
            out = 'z(t)'
            wen = '     y(t) = %s(z(t))' %(M.out.type)
            nlo = 1
        else:
            out = 'y(t)'
            wen = ''
            nlo = 0
        #endif
    else:
        out = 'y(t)'
        wen = ''
        nlo = 0
    #endif

    # Switch according to model type
    match M.type.lower():
        case 'nonpar':
            ep.alg = 'nonpar'
            ep.modelEquations = 'G(%s) = U(%s) / Y(%s)' %(M.op, M.op, M.op)


        case 'ar' | 'arx' | 'farx':
            match Z.type:
                case 'time':
                    ep.alg = 'barx'
                    for k in range(0, M.nu):
                        if not M.inp[k].type == 'linear':
                            ep.startNL = 'startNL'
                        #endif
                    #endif
                    for k in range(0, M.ny):
                        if not M.out.type == 'linear':
                            ep.startNL = 'startNL'
                        #endif
                    #endif
                case 'frequency':
                    ep.alg = 'farx'
            #endmatch

            if M.type == 'ar':
                cr = '\n'
                s1 = '                          '
                s2 = '     A(%s)y(t) = e(t)      ' %(M.op)
                s3 = '                          '
                s4 = 'Order of A: ' + str(M.nA[:])
                spc = '       '
                ep.modelEquations = cr + s1 + spc + s4 + cr + s2 + cr + s3 + cr
            elif M.type in ['arx', 'farx']:
                cr = '\n'
                s1 = '                                   '
                s2 = '     A(%s)y(t) = B(%s)u(t) + e(t)' %(M.op, M.op)
                s3 = '                                   '
                s4 = 'Order of A: ' + str(M.nA[:])
                s5 = 'Order of B: ' + str(M.nB[:])
                spc = '       '
                ep.modelEquations = cr + s1 + spc + s4 + cr + s2 + cr + s3 + spc + s5 + cr
            #endif


        case 'fir' | 'nfir':
            ep.alg = 'fir'
            for k in range(0, M.nu):
                if not M.inp[k].type == 'linear':
                    ep.startNL   = 'startNL'
                    ep.startG    = 'startG'
                    ep.alg       = OPT.alg
                #endif
            #endif
            for k in range(0, M.ny):
                if not M.out.type == 'linear':
                    ep.startNL   = 'startNL'
                    ep.startG    = 'startG'
                    ep.alg       = OPT.alg
                #endif
            #endif
            if M.nu > 1:
                g1 = ' %i               ' %(M.nu)
                g2 = 'sum  B_i(%s)%s' %(M.op, inp)
                g3 = 'i=1              '
            else:
                g1 = '          '
                g2 = 'B(%s)%s' %(M.op, inp)
                g3 = '          '
            #endif
            h1 =         '         '
            h2 = ' + e(t)'
            h3 =         '         '
            cr = '\n'
            s1 = '            ' + g1 + h1
            s2 = '      %s = ' %(out) + g2 + h2
            s3 = '            ' + g3 + h3
            s4 = 'Order of B: ' + str(M.nB[:].T)
            spc = '       '
            ep.modelEquations = cr + ham + cr + s1 + cr + s2 + spc + s4 + cr + s3 + cr + cr + wen + cr


        case 'static':
            ep.alg = OPT.alg
            cr = '\n'
            if M.nu > 1:
                g1 = ' %i               ' %(M.nu)
                g2 = 'sum %s' %(M.op, inp)
                g3 = 'i=1              '
            else:
                g1 = '    '
                g2 = '%s' %(inp)
                g3 = '    '
            #endif
            cr = '\n'
            s1 = '            ' + g1
            s2 = '     %s = ' %(out) + g2
            s3 = '            ' + g3
            ep.modelEquations = cr + ham + cr + s1 + cr + s2 + cr + s3 + cr
            for k in range(0, M.nu):
                if not M.inp[k].type == 'linear':
                    ep.startNL = 'startNL'
                #endif
            #endif
            for k in range (0, M.ny):
                if not M.out.type == 'linear':
                    ep.startNL = 'startNL'
                #endif
            #endif


        case 'arma' | 'armax' | 'oe' | 'bj':
            if not M.type == 'arma':
                ep.startG = 'startG'
            #endif
            ep.startH = 'startH'
            ep.alg = OPT.alg
            for k in range(0, M.nu):
                if not M.inp[k].type == 'linear':
                    ep.startNL = 'startNL'
                #endif
            #endfor
            if not M.out.type == 'linear':
                ep.startNL = 'startNL'
            #endif

            if M.type == 'arma':
                cr = '\n'
                s1 = '                          '
                s2 = '     A(%s)y(t) = C(%s)e(t)' %(M.op, M.op)
                s3 = '                          '
                s4 = 'Order of A: ' + str(M.nA)
                s5 = 'Order of C: ' + str(M.nC)
                spc = '       '
                ep.modelEquations = cr + s1 + spc + s4 + cr + s2 + cr + s3 + spc + s5 + cr
            elif M.type == 'armax':
                if M.nu > 1:
                    g1 = ' %i               ' %(M.nu)
                    g2 = 'sum  B_i(%s)%s' %(M.op,inp)
                    g3 = 'i=1              '
                else:
                    g1 = '          ' %(M.op)
                    g2 = 'B(%s)%s' %(M.op,inp)
                    g3 = '          ' %(M.op)
                #endif
                h1 =         '            '
                h2 = ' + C(%s)e(t)' %(M.op)
                h3 =         '            '
                cr = '\n'
                s1 = '            ' + g1 + h1
                s2 = ' A(%s)%s = ' %(M.op,out) + g2 + h2
                s3 = '            ' + g3 + h3
                s4 = 'Order of B: ' + str(M.nB[:].T) + '    Order of C: ' + str(M.nC)
                s5 = 'Order of A: ' + str(M.nA[:].T)
                spc = '       '
                ep.modelEquations = cr + ham + cr + s1 + spc + s4 + cr + s2 + cr + s3 + spc + s5 + cr + cr + wen + cr
            elif M.type == 'oe':
                if M.nu > 1:
                    g1 = ' %i   B_i(%s)       ' %(M.nu,M.op)
                    g2 = 'sum  ------%s' %(inp)
                    g3 = 'i=1  A_i(%s)       ' %(M.op)
                else:
                    g1 = 'B(%s)     ' %(M.op)
                    g2 = '----%s' %(inp)
                    g3 = 'A(%s)     ' %(M.op)
                #endif
                h1 = '       '
                h2 = ' + e(t)'
                h3 = '       '
                cr = '\n'
                s1 = '            ' + g1 + h1
                s2 = '     %s = ' %(out) + g2 + h2
                s3 = '            ' + g3 + h3
                s4 = 'Order of B: ' + str(M.nB[:].T)
                s5 = 'Order of A: ' + str(M.nA[:].T)
                spc = '       '
                ep.modelEquations = cr + ham + cr + s1 + spc + s4 + cr + s2 + cr + s3 + spc + s5 + cr + cr + wen + cr
            elif M.type == 'bj':
                if M.nu > 1:
                    g1 = ' %i   B_i(%s)       ' %(M.nu,M.op)
                    g2 = 'sum  ------%s' %(inp)
                    g3 = 'i=1  A_i(%s)       ' %(M.op)
                else:
                    g1 = 'B(%s)     ' %(M.op)
                    g2 = '----%s' %(inp)
                    g3 = 'A(%s)     ' %(M.op)
                #endif
                h1 = '   C(%s)     ' %(M.op)
                h2 = ' + ----e(t)'
                h3 = '   D(%s)     ' %(M.op)
                cr = '\n'
                s1 = '            ' + g1 + h1
                s2 = '     %s = ' %(out) + g2 + h2
                s3 = '            ' + g3 + h3
                s4 = 'Order of B: ' + str(M.nB[:].T) + '    Order of C: ' + str(M.nC)
                s5 = 'Order of A: ' + str(M.nA[:].T) + '    Order of D: ' + str(M.nD)
                spc = '       '
                ep.modelEquations = cr + ham + cr + s1 + spc + s4 + cr + s2 + cr + s3 + spc + s5 + cr + cr + wen + cr
            #endif


        case 'ss' | 'bilin' | 'bilinear' | 'lpv':
            match OPT.alg:
                case 'sid' | 'n4sid' | 'cca' | 'subspace':
                    ep.alg = 'subspace'

                case 'gn' | 'em':
                    ep.startG = 'startG'
                    ep.startH = 'startH'
                    ep.alg = OPT.alg
            #endmatch

            # Set the model equations
            if M.op == 's':
                s0 = '     u(s)   = u(t)  t <= s < t+d     assume piecewise constant input'
                s1 = '     .                         .   '
                s2 = '     x(t)   = Ax(t) + Bu(t) + Ke(t)'
                s3 = '     .                         .   '
                s4 = '     z(t)   = Cx(t)         +  e(t)'
                s5 = '               t+d .'
                s6 = '     y(t+d) =  int z(s)ds            assume integrated sampling '
                s7 = '                t '
                cr = '\n'
                ep.modelEquations = cr + s0 + cr + cr + s1 + cr + s2 + cr + cr + s3 + cr + s4 + cr + cr + s5 + cr + s6 + cr + s7 + cr + cr
            else:
                if M.nu > 0:
                    Bu = ' + Bu(t)'
                    if M.estD:
                        Du = ' + Du(t)'
                    else:
                        Du = ''
                    #endif
                    if M.estF:
                        Fukx = ' + F*kron(u(t),x(t))'
                    else:
                        Fukx = ''
                    #endif
                    if M.estG:
                        Gukx = ' + G*kron(u(t),x(t))'
                    else:
                        Gukx = ''
                    #endif
                else:
                    Fukx = ''
                    Gukx = ''
                    Bu   = ''
                    Du   = ''
                #endif
                if M.estK:
                    Ke = ' + Ke(t)'
                else:
                    Ke = ''
                #endif
                if M.type in ['bilin','bilinear']:
                    ep.modelEquations = \
                    '\n     %sx(t) = Ax(t)' %(M.op) + Bu + Fukx + Ke + '   ' + \
                    'nx = %i,  nu = %i,  ny = %i' %(M.nx, M.nu, M.ny) + \
                    '\n      y(t) = Cx(t)' + Du + Gukx + \
                    ' +  e(t)\n' %(M.op)
                else:
                    ep.modelEquations = \
                    '\n     %sx(t) = Ax(t)' %(M.op) + Bu + Ke + '   ' + \
                    'nx = %i,  nu = %i,  ny = %i' %(M.nx, M.nu, M.ny) + \
                    '\n      y(t) = Cx(t)' + Du + \
                    ' +  e(t)\n'
                #endif
            #endif


        case 'nlss':
            # Set empty model equations
            ep.modelEquations = ''

            # Make est call the main routine emnlss
            ep.alg     = 'emnlss'
            ep.finishM = 'nlssfinish'


        case _:
            raise Exception("M.type is unknown!")

    #endmatch

    if 'modelEquations' not in ep:
        ep.modelEquations = ''
    #endif


    return ep

