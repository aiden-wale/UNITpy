
import numpy as np
import uonidtoolbox as unit

length = unit._utils.length
isempty = unit._utils.isempty


def startM(*args):

    nargin = len(args)

    # ============================== Get Z and M ===============================
    match nargin:
        case 0:
            Z = unit.struct()
            M = unit.struct()
            nu = 1
            ny = 1

        case 1: # input should be pertaining to M
            Z = unit.struct()
            M = unit.struct()
            nu = 1
            ny = 1

            if isinstance(args[0], str):
                M.type = args[0]
            elif isinstance(args[0], (unit.struct, dict)):
                if 'type' in args[0]:
                    if args[0]['type'] in ['time', 'frequency']:
                        # then args[0] pertains to Z
                        Z = unit.struct(args[0])
                        M = unit.struct()
                    else:
                        M = unit.struct(args[0])
                        if 'nu' in M: nu = M.nu
                        if 'ny' in M: ny = M.ny
                    #endif
                else:
                    M = args[0]
                    if 'nu' in M: nu = M.nu
                    if 'ny' in M: ny = M.ny
                #endif
            elif isinstance(args[0], (int, float, complex)):
                M.nA = args[0]
            #endif

        case _: # 2 args
            Z = unit.startZ(args[0])
            nu = Z.nu
            ny = Z.ny

            if isinstance(args[1], unit.struct):
                M = args[1]
            if isinstance(args[1], dict):
                M = unit.struct(args[1])
            elif isinstance(args[1], str):
                M = unit.struct()
                M.type = args[1]
            elif isinstance(args[1], (int, float, complex)):
                M = unit.struct()
                if args[1] >= 0:
                    M.nA = args[1]
                #endif
            #endif
            if 'nu' in M: nu = M.nu
            if 'ny' in M: ny = M.ny
    #endmatch
    # ==========================================================================

    M = _inputCleanse_M(M)
    if 'nu' in M: nu = M.nu
    if 'ny' in M: ny = M.ny


    # ================== THE DEFAULT MODEL STRUCTURE ===========================
    m = unit.struct()
    gord        = 5         # Default order of G dynamics
    hord        = 2         # Default order of H dynamics
    m.type   = 'arx'     # Default type of model

    # Transfer function defaults
    m.A      = np.array([]).reshape([0,0])
    m.B      = np.array([]).reshape([0,0])
    m.C      = np.array([]).reshape([0,0])
    m.D      = np.array([]).reshape([0,0])
    m.nA     = gord*np.ones([nu])
    m.nB     = gord*np.ones([nu])
    m.nC     = hord*np.ones([1])
    m.nD     = hord*np.ones([1])

    # State-space defaults
    m.nx     = gord*np.ones([1])
    m.estD   = 1
    m.estK   = 1
    m.estX1  = 1
    m.estF   = 1
    m.estG   = 1
    m.par    = 'ddlc'

    # Other deafults
    m.theta      = np.array([]).reshape([0,0])
    m.nu         = nu
    m.ny         = ny
    m.finishM    = 'finishM'

    # If a sample time was given by the data, then use it
    if 'T' in Z:
        m.T = Z.T
    else:
        m.T = 1
    #endif

    # Default operator depends on the type of data (time domain = 'q', freq. dom. = 's')
    if 'type' in Z:
        match Z.type:
            case 'frequency':
                m.op = 's'
            case _:
                m.op = 'q'
        #endmatch
    else:
        m.op = 'q'
    #endif

    if 'w' in Z:
        m.w = Z.w.reshape(1, Z.w.size)
    else:
        m.w = np.logspace(np.log10(np.pi/m.T/1000), np.log10(np.pi/m.T), 1000).reshape(1, 1000)
    #endif
    m.delay = np.zeros([np.max([nu,1])])
    # ==========================================================================


    # ===================== WHAT IS THE MODEL TYPE? ============================
    # If type not specified then guess it
    if 'type' not in M:
        
        # Determine what the user has supplied in terms of poly's, orders and state-space matrices
        isApoly = 'A' in M or 'nA' in M
        isBpoly = 'B' in M or 'nB' in M
        isCpoly = 'C' in M or 'nC' in M
        isDpoly = 'D' in M or 'nD' in M
        if 'ss' in M:
            isAss = 'A' in M.ss
            isBss = 'B' in M.ss
            isCss = 'C' in M.ss
        else:
            isAss = 0
            isBss = 0
            isCss = 0
        #endif

        # Based on supplied orders and/or poly's, then try and guess type
        if nu == 0:
            if      isApoly and  not isBpoly and  not isCpoly and  not isDpoly: M.type = 'ar'
            if      isApoly and  not isBpoly and      isCpoly and  not isDpoly: M.type = 'arma'
        #endif
        if  not isApoly and      isBpoly and  not isCpoly and  not isDpoly: M.type = 'fir'
        if      isApoly and      isBpoly and  not isCpoly and  not isDpoly: M.type = 'arx'
        if      isApoly and      isBpoly and      isCpoly and  not isDpoly: M.type = 'armax'
        if      isApoly and      isBpoly and  not isCpoly and  not isDpoly: M.type = 'oe'
        if      isApoly and      isBpoly and      isCpoly and      isDpoly: M.type = 'bj'
        # If a D order or poly is given, then must be BJ form
        if                                                         isDpoly: M.type = 'bj'
        # If missing any poly information, then look to a SS type
        if  not isApoly and  not isBpoly and  not isCpoly and  not isDpoly:
            if  isAss or isBss or isCss:
                M.type = 'ss'
            #endif
        #endif

        # If there are multiple outputs, then default ot a state-space model
        if ny > 1:
            M.type = 'ss'
        #endif

        # Fill in default type if not yet specified
        if 'type' not in M:
            M.type = m.type
        #endif
    #endif
    # ==========================================================================


    # ========================= FILL IN DEFAULTS =============================== 

    # ------------------------------- AR ---------------------------------------
    if M.type in ['ar', 'nar']:
        # Has the user specified an order for the A polynomial?
        if 'nA' not in M:
            if 'A' not in M:
                M.nA = gord*np.ones([ny])
            else:
                M.nA = np.array([length(M.A)-1])
            #endif
        #endif

        # Make sure that the nA and A variables have only one entry
        # M.nA = M.nA[0]
        if length(M.nA) > 1: M.nA = M.nA[0]
        if 'A' in M:
            if not isempty(M.A):
                M.A = np.array(M.A)
                M.A = M.A[0,:].reshape([1,M.A.shape[1]])
            #endif
        #endif
        M.B = np.array([[0.0]])
        M.C = np.array([[1.0]])
        M.D = np.array([[1.0]])

        # Set default orders for B, C and D poly's
        M.nB = np.array([0])
        M.nC = np.array([0])
        M.nD = np.array([0])

        # Make sure the number of inputs == 0
        M.nu = 0
        M.ny = 1
    #endif
    # --------------------------------------------------------------------------


    # ------------------------------ ARMA --------------------------------------
    if M.type in ['arma', 'narma']:
        # Has the user specified an order for the A polynomial?
        if 'nA' not in M:
            if 'A' not in M:
                M.nA = gord
            else:
                if np.floor(M.A) == M.A:
                    M.nA = M.A
                else:
                    M.nA = length(M.A)-1
                #endif
            #endif
        #endif
        # Has the user specified an order for the C polynomial?
        if 'nC' not in M:
            if 'C' not in M:
                M.nC = M.nA
            else:
                if np.floor(M.C) == M.C or length((M.C)) < 2:
                    M.nC = M.C
                else:
                    M.nC = length(M.C)-1
                #endif
            #endif
        #endif

        # Make sure that the nA, nC and A,C variables have only one entry
        if length(M.nA) > 1: M.nA = M.nA[0]
        if 'A' in M:
            if not isempty(M.A):
                M.A = np.array(M.A)
                M.A = M.A[0,:].reshape([1,M.A.shape[1]])
            #endif
        #endif
        if length(M.nC) > 1: M.nC = M.nC[0]
        if 'C' in M:
            if not isempty(M.C):
                M.C = np.array(M.C)
                M.C = M.C[0,:].reshape([1,M.C.shape[1]])
            #endif
        #endif
        M.B = np.array([[0.0]])
        M.D = np.array([[1.0]])

        # Set default orders for B, C and D poly's
        M.nB = np.array([0])
        M.nD = np.array([0])

        # Make sure the number of inputs == 0
        M.nu = 0
        M.ny = 1
    #endif
    # --------------------------------------------------------------------------


    # ------------------------------ FIR ---------------------------------------
    if M.type in ['fir', 'nfir']:
        # This may also include generalised FIR models for onid
        # Has the user specified an order for the B polynomial?
        if 'nB' not in M:
            if 'B' not in M:
                M.nB = gord*np.ones([nu,1])
                M.nA = np.ones([nu,1])
            else:
                for i in range(0,nu):
                    M.nA = np.zeros([nu,1])
                    if np.floor(M.B) == M.B:
                        M.nB[i] = M.B[i]
                    else:
                        M.nB[i] = length(M.B[i,:])-1
                    #endif
                #endfor
            #endif
        #endif
        #If length(M.nB) != number of inputs, then extend it so that it does
        if length(M.nB) < nu:
            df = nu - length(M.nB)
            # M.nB = np.vstack((M.nB, M.nB[-1]*np.ones([df,1])))
            if length(M.nB) > 1:
                M.nB = np.vstack((M.nB, M.nB[-1]*np.ones([df,1])))
            else:
                M.nB = np.vstack((M.nB, M.nB*np.ones([df,1])))
            #endif
        else:
            M.nu = length(M.nB)
        #endif
        #If length(M.nA) != number of inputs, then extend it so that it does
        M.nA = 0
        if length(M.nA) != nu:
            df = nu - length(M.nA)
            if df > 0:
                # M.nA = np.vstack((M.nA, M.nA[-1]*np.ones([df,1])))
                if length(M.nA) > 1:
                    M.nA = np.vstack((M.nA, M.nA[-1]*np.ones([df,1])))
                else:
                    M.nA = np.vstack((M.nA, M.nA*np.ones([df,1])))
                #endif
            elif df < 0:
                M.nA = M.nA[0:-df]
            #endif
        #endif
        M.nC = 0
        M.nD = 0
        if 'poles' not in M:
            M.poles = np.zeros([1, int(np.max(M.nB))]) # Special case of FIR from onid.m
        #endif

        # Set all remaining polynomials to default values
        if 'A' not in M: M.A = np.ones([nu,1])
        M.C = np.array([[1.0]])
        M.D = np.array([[1.0]])
    #endif
    # --------------------------------------------------------------------------


    # ------------------------------ ARX ---------------------------------------
    if M.type in ['arx', 'narx']:
        # Has the user specified an order for the A polynomial?
        if 'nA' not in M:
            if 'A' not in M:
                M.nA = gord
            elif isempty(M.A):
                M.nA = gord
            else:
                if M.A.size == 1:
                    if np.floor(M.A) == M.A:
                        M.nA = M.A
                    #endif
                else:
                    M.nA = length(M.A)-1
                #endif
            #endif
        #endif
        # Has the user specified an order for the B polynomial?
        if 'nB' not in M:
            if 'B' not in M:
                M.nB = M.nA
            elif isempty(M.B):
                M.nB = M.nA
            else:
                if M.B.size == 1:
                    if np.floor(M.B) == M.B:
                        M.nB = M.B
                    #endif
                else:
                    for i in range(0,nu):
                        M.nB = length(M.B)-1
                    #endfor
                #endif
            #endif
        #endif
        # If length(M.nB)~=number of inputs, then extend it so that it does
        if length(M.nB) < nu:
            df = nu - length(M.nB)
            # M.nB = np.vstack([M.nB, M.nB[-1]*np.ones([df,1])])
            if length(M.nB) > 1:
                M.nB = np.vstack((M.nB, M.nB[-1]*np.ones([df,1])))
            else:
                M.nB = np.vstack((M.nB, M.nB*np.ones([df,1])))
            #endif
        else:
            M.nu = length(M.nB)
        #endif

        # Make sure that the nA and A variables have only one entry
        # M.nA = M.nA[0]
        if length(M.nA) > 1: M.nA = M.nA[0]
        if 'A' in M:
            if not isempty(M.A):
                # M.A = np.array(M.A)
                M.A = M.A[0,:]
            #endif
        #endif
        M.C = np.array([[1.0]])
        M.D = np.array([[1.0]])

        # Set default orders for C and D poly's
        M.nC = np.array([0])
        M.nD = np.array([0])
    #endif
    # --------------------------------------------------------------------------


    # ----------------------------- ARMAX --------------------------------------
    if M.type in ['armax', 'narmax']:
        # Has the user specified an order for the A polynomial?
        if 'nA' not in M:
            if 'A' not in M:
                M.nA = gord
            else:
                if M.A.size == 1:
                    if np.floor(M.A) == M.A:
                        M.nA = M.A
                    #endif
                else:
                    M.nA = length(M.A)-1
                #endif
            #endif
        #endif
        # Has the user specified an order for the B polynomial?
        if 'nB' not in M:
            if 'B' not in M:
                M.nB = M.nA
            else:
                M.nB = np.zeros(nu)
                for i in range(0,nu):
                    if np.floor(M.B) == M.B:
                        M.nB[i] = M.B[i]
                    else:
                        M.nB[i] = length(M.B)-1
                    #endif
                #endfor
            #endif
        #endif
        # Has the user specified an order for the C polynomial?
        if 'nC' not in M:
            if 'C' not in M:
                M.nC = hord
            else:
                if np.floor(M.C) == M.C:
                    M.nC = M.C
                else:
                    M.nC = length(M.C)-1
                #endif
            #endif
        #endif
        # If length(M.nB)~=number of inputs, then extend it so that it does
        if length(M.nB) < nu:
            df = nu - length(M.nB)
            # M.nB = np.vstack((M.nB, M.nB[-1]*np.ones([df,1])))
            if length(M.nB) > 1:
                M.nB = np.vstack((M.nB, M.nB[-1]*np.ones([df,1])))
            else:
                M.nB = np.vstack((M.nB, M.nB*np.ones([df,1])))
            #endif
        else:
            M.nu = length(M.nB)
        #endif

        # Make sure that the nA and A variables have only one entry
        # M.nA = M.nA[0] 
        if length(M.nA) > 1: M.nA = M.nA[0]
        if 'A' in M:
            if not isempty(M.A):
                M.A = np.array(M.A)
                M.A = M.A[0,:]
            #endif
        #endif
        M.D = np.array([[1.0]])

        # Set default orders for C and D poly's
        M.nD = np.array([0])
    #endif
    # --------------------------------------------------------------------------


    # ------------------------------ OE ----------------------------------------
    if M.type in ['oe', 'noe']:
        # Has the user specified an order for the A polynomial?
        if 'nA' not in M:
            if 'A' not in M:
                M.nA = gord
            else:
                M.nA = np.zeros(nu)
                for i in range(0,nu):
                    if np.floor(M.A) == M.A:
                        M.nA[i] = M.A[i]
                    else:
                        M.nA[i] = length(M.A[i,:])-1
                    #endif
                #endfor
            #endif
        #endif
        # Has the user specified an order for the B polynomial?
        if 'nB' not in M:
            if 'B' not in M:
                M.nB = M.nA
            else:
                M.nB = np.zeros(nu)
                for i in range(0,nu):
                    if np.floor(M.B) == M.B:
                        M.nB[i] = M.B[i]
                    else:
                        M.nB[i] = length(M.B[i,:])-1
                    #endif
                #endfor
            #endif
        #endif
        # If length(M.nA)~=number of inputs, then extend it so that it does
        if length(M.nA) < nu:
            df = nu - length(M.nA)
            # M.nA = np.vstack((M.nA, M.nA[-1]*np.ones([df,1])))
            if length(M.nA) > 1:
                M.nA = np.vstack((M.nA, M.nA[-1]*np.ones([df,1])))
            else:
                M.nA = np.vstack((M.nA, M.nA*np.ones([df,1])))
            #endif
        else:
            M.nu = length(M.nA)
        #endif
        # If length(M.nB)~=number of inputs, then extend it so that it does
        if length(M.nB) < nu:
            df = nu - length(M.nB)
            # M.nB = np.vstack((M.nB, M.nB[-1]*np.ones([df,1])))
            if length(M.nB) > 1:
                M.nB = np.vstack((M.nB, M.nB[-1]*np.ones([df,1])))
            else:
                M.nB = np.vstack((M.nB, M.nB*np.ones([df,1])))
            #endif
        else:
            M.nu = length(M.nB)
        #endif

        # Set default orders for C and D poly's
        M.nC = np.array([0])
        M.nD = np.array([0])
        M.C  = np.array([[1.0]])
        M.D  = np.array([[1.0]])
    #endif
    # --------------------------------------------------------------------------


    # ------------------------------ BJ ----------------------------------------
    if M.type in ['bj', 'nbj']:
        # Has the user specified an order for the A polynomial?
        if 'nA' not in M:
            if 'A' not in M:
                M.nA = gord
            else:
                M.nA = np.zeros(M.A.shape[0])
                for i in range(0, M.A.shape[0]):
                    if np.floor(M.A) == M.A:
                        M.nA[i] = M.A[i]
                    else:
                        M.nA[i] = length(M.A[i,:])-1
                    #endif
                #endfor
            #endif
        #endif
        # Has the user specified an order for the B polynomial?
        if 'nB' not in M:
            if 'B' not in M:
                M.nB = M.nA
            else:
                M.nB = np.zeros(M.B.shape[0])
                for i in range(0, M.B.shape[0]):
                    if np.floor(M.B) == M.B:
                        M.nB[i] = M.B[i]
                    else:
                        M.nB[i] = length(M.B[i,:])-1
                    #endif
                #endfor
            #endif
        #endif
        # Has the user specified an order for the C polynomial?
        if 'nC' not in M:
            if 'C' not in M:
                if 'D' in M:
                    if np.floor(M.D) == M.D: hord = M.D
                #endif
                M.nC = hord
            else:
                if np.floor(M.C) == M.C and length(M.C) < 2:
                    M.nC = M.C
                else:
                    M.nC = length(M.C)-1
                #endif
            #endif
        #endif
        # Has the user specified an order for the D polynomial?
        if 'nD' not in M:
            if 'D' not in M:
                M.nD = M.nC
            else:
                if np.floor(M.D) == M.D:
                    M.nD = M.D
                else:
                    M.nD = length(M.D)-1
                #endif
            #endif
        #endif
        # If length(M.nA)~=number of inputs, then extend it so that it does
        if length(M.nA) < nu:
            df = nu - length(M.nA)
            # M.nA = np.vstack((M.nA, M.nA[-1]*np.ones([df,1])))
            if length(M.nB) > 1:
                M.nA = np.vstack((M.nA, M.nA[-1]*np.ones([df,1])))
            else:
                M.nA = np.vstack((M.nA, M.nA*np.ones([df,1])))
            #endif
        else:
            M.nu = length(M.nA)
        #endif
        # If length(M.nB)~=number of inputs, then extend it so that it does
        if length(M.nB) < nu:
            df = nu - length(M.nB)
            # M.nB = np.vstack((M.nB, M.nB[-1]*np.ones([df,1])))
            if length(M.nB) > 1:
                M.nB = np.vstack((M.nB, M.nB[-1]*np.ones([df,1])))
            else:
                M.nB = np.vstack((M.nB, M.nB*np.ones([df,1])))
            #endif
        else:
            M.nu = length(M.nB)
        #endif
    #endif
    # --------------------------------------------------------------------------


    # ------------------------------ SS ----------------------------------------
    if M.type in ['ss', 'nss', 'bilin', 'bilinear']:
        # Determine state order nx
        if 'nx' not in M:
            if 'ss' in M:
                if 'A' in M.ss:
                    M.nx = M.ss.A.shape[0]
                #endif
            #endif
        #endif

        if 'nx' not in M:
            if 'nA' in M:
                M.nx = M.nA
            elif 'A' in M:
                if np.floor(M.A) == M.A:
                    M.nx = np.max(M.A) # TODO: Check this (supremum of matrix elements?)
                else:
                    M.nx = m.nx
                #endif
            else:
                M.nx = m.nx
            #endif
        #endif

        # If a system has been given, then make sure nx, nu, and ny are set correctly
        if 'ss' in M:
            if 'A' in M.ss:
                if not isempty(M.ss.A):
                    M.nx = M.ss.A.shape[0]
                #endif
            #endif
            if 'B' in M.ss:
                if not isempty(M.ss.B):
                    if nu != M.ss.B.shape[1]:
                        M.nu = M.ss.B.shape[1]
                    #endif
                #endif
            else:
                M.nu = 0
            #endif
            if 'C' in M.ss:
                if not isempty(M.ss.C):
                    if ny != M.ss.C.shape[0]:
                        M.ny = M.ss.C.shape[0]
                    #endif
                #endif
            else:
                M.ny = 0
            #endif
        #endif

        # Make sure that if a continuous-time model is asked for then we set M.par = 'full'
        discrete = 1
        if 'op' in M:
            if M.op == 's':
                discrete = 0
            #endif
        #endif
        if 'par' not in M:
            if not discrete:
                M.par = 'full'
            #endif
        #endif
    #endif
    # --------------------------------------------------------------------------


    # ---------------------------- NONPAR --------------------------------------
    if M.type in ['nonpar']:
        pass
    #endif
    # --------------------------------------------------------------------------


    # ---------------------------- STATIC --------------------------------------
    if M.type in ['static']:
        M.A  = np.array([[0.0]])
        M.B  = np.array([[0.0]])
        M.C  = np.array([[0.0]])
        M.D  = np.array([[0.0]])
        M.nA = np.array([1])
        M.nB = np.array([0])
        M.nC = np.array([1])
        M.nD = np.array([1])
    #endif
    # --------------------------------------------------------------------------


    # -------------------- NONLINEAR STATE-SPACE -------------------------------
    if M.type in ['nlss']:
        M.op         = 'q'
        M.finishM    = 'nlssfinish'
    #endif
    # --------------------------------------------------------------------------
    

    # --------------------------- DEFAULT --------------------------------------
    # case _:
    #     M.type = m.type

    # ------------------Fill in bits of M that come from data-------------------
    if 'type' in Z:
        match Z.type:
            case 'frequency':
                if 'w' not in M:
                    M.w = Z.w
                #endif
                if 'T' not in M:
                    if M.T*np.max(M.w) > np.pi + 1000*unit._utils.eps:
                        # warning('M.T and Z.w are not compatible for M.op="q", resetting M.T = pi/max(M.w);')
                        M.T = np.pi/np.max(M.w)
                    #endif
                else:
                    M.T = np.pi/np.max(M.w)
                #endif

            case 'time':
                if 'T' not in M:
                    M.T = Z.T
                #endif
                if 'w' not in M:
                    if 'w' in Z:
                        M.w = Z.w
                    else:
                        M.w = np.logspace(np.log10(np.pi/M.T)-4, np.log10(np.pi/M.T), 2000).reshape(1, 2000)
                    #endif
                #endif

            case _:
                # shouldnt get here
                pass

            #endmatch
        #endif
    #endif
    # --------------------------------------------------------------------------


    # ------------------Fill in missing or empty fields of M--------------------
    if not M.type in ['nlss']: # If M.type is nonlinear state-space then don't fill in defaults
        if not M:
            M = m.copy()
        else:
            for k in m.keys():
                if k not in M:
                    M[k] = m[k]
                elif isempty(M[k]):
                    M[k] = m[k]
                #endif
            #endif
        #endif

        # Make sure the delays are of the correct orientation and size.
        M.delay = M.delay.reshape(M.delay.size)
        if M.delay.size != M.nu and M.nu > 0:
            # M.delay = np.hstack([M.delay, np.zeros([np.max([0, M.nu - M.delay.shape[0]]), 1])])
            M.delay = np.array([M.delay, np.zeros([M.nu-M.delay.size])]).reshape(nu)
        #endif
        M.delay = np.array(M.delay, dtype='int')

        # Initialise Hammerstein and Wiener NL blocks
        M = unit.startNL(Z,M)

    #endif

    M = unit.struct(dict(sorted(M.items())))

    return M;


def _inputCleanse_M(M):
    # Ensure any given M.{A,B,C,D} are 2D numpy arrays
    for k in ['A', 'B', 'C', 'D']:
        if k in M:
            if not isinstance(M[k], np.ndarray):
                if isempty(M[k]):
                    M[k] = np.array([[]])
                else:
                    M[k] = np.array([M[k]]).reshape(1,1)
                #endif
            elif M[k].ndim == 1:
                M[k] = M[k].reshape(1, M[k].size)
            elif M[k].ndim > 2:
                if np.prod(M[k].shape) != np.max(M[k].shape):
                    raise Exception(f"M.{k} must be a 2D numpy array")
                #endif
            #endif
        #endif
    #endfor

    # Ensure any given M.{nA,nB,nC,nD} are 1D numpy arrays
    for k in ['nA', 'nB', 'nC', 'nD']:
        if k in M:
            if not isinstance(M[k], np.ndarray):
                if isempty(M[k]):
                    M[k] = np.array([])
                else:
                    M[k] = np.array([M[k]]).reshape(1,1) # TODO: these are meant to be 1D...
                #endif
            elif M[k].ndim > 1:
                if np.prod(M[k].shape) != np.max(M[k].shape):
                    raise Exception(f"M.{k} must be a 1D numpy array")
                #endif
            #endif
            M[k] = np.array(M[k], dtype='int')
        #endif
    #endfor

    if 'w' in M:
        M.w = np.squeeze(M.w)
        if len(M.w.shape) > 1: raise Exception("M.w should be a 1-D array")
    #endif
    if 'delay' in M:
        M.delay = np.array([M.delay], dtype='int')
        M.delay = M.delay.reshape(M.delay.size)
    #endif

    for k in ['estD','estF','estG','estK','estX1','nu','nx','ny']:
        if k in M:
            if isinstance(M[k], np.ndarray):
                if M[k].size != 1: 
                    raise Exception(f"Model variable: M.{k} must be a scalar integer")
                #endif
                M[k] = int(M[k].flatten()[0])
            else:
                M[k] = int(M[k])
            #endif
        #endif
    #endfor

    return M
#endfunction


def _make2d_SS_matrices(SS,nx,nu,ny):

    for k in ['A','B','C','D','G','R','S']:
        if not isinstance(SS[k], np.ndarray):
            SS[k] = np.array([SS[k]])
        #endif
    #endfor

    if SS.A.size > 0: SS.A = SS.A.reshape(nx,nx)
    if SS.B.size > 0: SS.B = SS.B.reshape(nx,nu)
    if SS.C.size > 0: SS.C = SS.C.reshape(ny,nx)
    if SS.D.size > 0: SS.D = SS.D.reshape(ny,nu)
    if SS.G.size > 0: SS.G = SS.G.reshape(nx,nx)
    if SS.R.size > 0: SS.R = SS.R.reshape(ny,ny)
    if SS.S.size > 0: SS.S = SS.S.reshape(nx,ny)

    return SS
#endfunction

