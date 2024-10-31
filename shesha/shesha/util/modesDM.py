import numpy as np



def computePistonFilteringMatrix(n):
    '''
    Creates a sqaure matrix that removes (=filter out) piston mode.
    <n> : size of the matrix

    P = computePistonFilteringMatrix(64)
    '''
    P = np.full((n, n), -1. / n)
    for i in range(n):
        P[i, i] += 1.0
    return P


def applySymcode2Coord(x, y, symcode=0):
    """
    Change signs of x or y, and/or swap them, according to symcode (0 to 7)
    <x>, <y>  : ndarrays. Coordinates.
    <symcode> : int, optional. Symcode, 0 to 7. The default is 0.

    Returns
    -------
    x, y : Transformed coordinates.
    """
    # Application of symcode
    symcode = symcode % 8
    symcode = int(symcode)
    if symcode>3:
        symcode = symcode - 4
        x, y = y, x
    if symcode==0:
        pass
    elif symcode==1:
        x *= -1
    elif symcode==2:
        y *= -1
    elif symcode==3:
        x *= -1
        y *= -1
    else:
        raise ValueError(f'Symcode value is {symcode}')
    return (x, y)



def generateActuCoordinates(nx, ny, Nactu, symcode=0):
    """
    Generates the coordinates of the actuators of the DM,
    centred on (0,0) with a pitch of 1.00.
    <nx>     : int. Number of actus across pupil diameter
    <ny>     : int. Number of actus across pupil diameter
    <Nactu>  : int. Total valid number of actus
    <symcode>: symcode, 0 to 7 (flip/mirror coords)

    nx, ny = 64, 64
    Nactu = 3228
    x, y, mask = generateActuCoordinates(nx, ny, Nactu)
    """
    # generated indices with pitch = 1.0
    x, y = np.indices((nx, ny))  # convention [x, y]
    # Centre things around 0
    x = x - np.mean(x.flatten())
    y = y - np.mean(y.flatten())
    # compute distance to centre
    r = x * x + y * y
    rf = r.flatten()
    # sort actuators by distances
    indr = np.argsort(rf)  # sort in increasing order
    rmax = (rf[indr[Nactu]] + rf[indr[Nactu - 1]]) / 2.
    x = x[r < rmax]
    y = y[r < rmax]
    # The symcode is applied to take into account the way the DM is seen on
    # the WFS. To be really clean, the same symcode operation should be applied
    # to <r> as well. However, <r> is so symmetric that the symcode transform
    # has no effect on it: the operation is skipped.
    x, y = applySymcode2Coord(x, y, symcode=symcode)
    return x, y, r < rmax


def geometricModes(x, y, coupling, filterPiston):
    '''
    Return modes
    <x, y>     : come from the output of the function
                 x, y, mask = generateActuCoordinates(nx, ny, Nactu)
    <coupling> : coupling factor between actus
    <filterPiston> : True/False, to filter the piston out or not

    nx, ny = 64, 64
    Nactu = 3228
    x, y, mask = generateActuCoordinates(nx, ny, Nactu)

    coupling = 0.4
    modes =  geometricModes(x, y, coupling, True)
    '''
    # Computation of the matrix of distances (squared !)
    print('Computing matrix of distances')
    mdist2 = (x[:, None] - x[None, :])**2 + (y[:, None] - y[None, :])**2
    print('Computing delta matrix')
    # delta = np.exp(mdist2 * np.log(coupling) )
    delta = 1. / (mdist2 * (1. / coupling - 1.) + 1.)
    if filterPiston:
        print('Filtering piston out')
        Nactu = delta.shape[0]
        P = computePistonFilteringMatrix(Nactu)
        delta = P.dot(delta.dot(P))
    print('Diagonalisation')
    l, U = np.linalg.eigh(delta)
    print('done')
    return U[:, ::-1]


def KLmodes(x, y, L0, filterPiston):
    '''
    Return modes
    <x, y>     : x,y coordinates of the DM actuators. 
    Could come from the output of the function:
    x, y, mask = generateActuCoordinates(nx, ny, Nactu) 
    OR can be any x,y DM coordinates. 

    <L0>       : outer scale, in the same units as x and y.
    <filterPiston> : True/False, to filter the piston out or not


    example of use:
    L0 = 25

    nx, ny = 64, 64
    Nactu = 3228
    x, y, mask = generateActuCoordinates(nx, ny, Nactu)

    modes =  KLmodes(x, y, L0, True)
    '''
    Nactu = len(x)
    # Computation of the matrix of distances (squared !)
    print('Computing matrix of distances')
    mdist2 = (x[:, None] - x[None, :])**2 + (y[:, None] - y[None, :])**2
    # normalisation of distances with respect to pupil diameter
    D = np.ptp(x)
    mdist2 /= D**2
    print('Computing covariance matrix')
    if L0 == 0:
        kolmo = -0.5 * 6.88 * mdist2**(5. / 6)
    else:
        kolmo = -0.5 * rodconan(np.sqrt(mdist2), L0 / D)
    if filterPiston:
        print('Filtering piston out')
        P = computePistonFilteringMatrix(Nactu)
        kolmo = P.dot(kolmo.dot(P))
    print('Diagonalisation')
    # L'unite des valuers propres est en radians^2 a la longueur d'onde ou
    # est exprime r0, et sachant que tout est normalise pour (D/r0)=1.
    #
    l, U = np.linalg.eigh(kolmo)
    l = l / Nactu
    print('done')
    return U[:, ::-1], l[::-1]

def macdo_x56(x,k=10):
    """   
    Computation of the Mc Donald function.

    f(x) = x**(5/6)*K_{5/6}(x)
    using a series for the esimation of K_{5/6}, taken from Rod Conan thesis :
    K_a(x)=1/2 \sum_{n=0}**\infty \frac{(-1)**n}{n!}
    \left(\Gamma(-n-a) (x/2)**{2n+a} + \Gamma(-n+a) (x/2)**{2n-a} \right) ,
    with a = 5/6.

    Setting x22 = (x/2)**2, setting uda = (1/2)**a, and multiplying by x**a,
    this becomes :
    x**a * Ka(x) = 0.5 $ -1**n / n! [ G(-n-a).uda x22**(n+a) + G(-n+a)/uda x22**n ]
    Then we use the following recurrence formulae on the following quantities :
    G(-(n+1)-a) = G(-n-a) / -a-n-1
    G(-(n+1)+a) = G(-n+a) /  a-n-1
    (n+1)! = n! * (n+1)
    x22**(n+1) = x22**n * x22
    and at each iteration on n, one will use the values already computed
    at step (n-1).
    The values of G(a) and G(-a) are hardcoded instead of being computed.

    The first term of the series has also been skipped, as it
    vanishes with another term in the expression of Dphi.
    
    """
    x = np.array(x) # Safe check
    a = 5./6.
    fn = 1.                             # initialisation factorielle 0!=1
    x2a = x**(2.*a)
    x22 = x*x/4.                        #  (x/2)**2
    x2n = 0.5                           # init (1/2) * x**0
    Ga  =  2.01126983599717856777       # Gamma(a) / (1/2)**a
    Gma = -3.74878707653729348337       # Gamma(-a) * (1/2.)**a
    s = np.zeros(x.shape)
    for n in range(k+1):
      dd = Gma * x2a
      if n:
        dd += Ga
      dd *= x2n
      dd /= fn
      # addition to s, with multiplication by (-1)**n
      if n%2:
          s -= dd
      else:
          s += dd
      # prepare recurrence iteration for next step
      if n<k:
        fn *= n+1     # factorial
        Gma /= -a-n-1 # gamma function
        Ga /= a-n-1   # idem
        x2n *= x22    # x**n
    return s



def asymp_macdo(x):
    """
    Computes a term involved in the computation of the phase struct
    function with a finite outer scale according to the Von-Karman
    model. The term involves the MacDonald function (modified bessel
    function of second kind) K_{5/6}(x), and the algorithm uses the
    asymptotic form for x ~ infinity.
    Warnings :
        - This function makes a floating point interrupt for x=0
    and should not be used in this case.
        - Works only for x>0.
    
    """
    x = np.array(x)
    # k2 is the value for
    # gamma_R(5./6)*2**(-1./6)
    k2 = 1.00563491799858928388289314170833
    k3 = 1.25331413731550012081   #  sqrt(pi/2)
    a1 = 0.22222222222222222222   #  2/9
    a2 = -0.08641975308641974829  #  -7/89
    a3 = 0.08001828989483310284   # 175/2187
    x_1 = 1./x
    res = k2 - k3*np.exp(-x)*x**(1/3.)*(1.0 + x_1*(a1 + x_1*(a2 + x_1*a3)))
    return res




def rodconan(r,L0,k=10):
    """ DOCUMENT rodconan(r,L0,k=)
    The phase structure function is computed from the expression
    Dphi(r) = k1  * L0**(5./3) * (k2 - (2.pi.r/L0)**5/6 K_{5/6}(2.pi.r/L0))

    For small r, the expression is computed from a development of
    K_5/6 near 0. The value of k2 is not used, as this same value
    appears in the series and cancels with k2.
    For large r, the expression is taken from an asymptotic form.
    
    """
    # k1 is the value of :
    # 2*gamma_R(11./6)*2**(-5./6)*pi**(-8./3)*(24*gamma_R(6./5)/5.)**(5./6)
    k1 = 0.1716613621245709486
    dprf0 = (2*np.pi/L0)*r
    # k2 is the value for gamma_R(5./6)*2**(-1./6),
    # but is now unused
    # k2 = 1.0056349179985892838    
    res = np.zeros(r.shape)
    Xlim = 0.75*2*np.pi
    largeX = dprf0>Xlim

    res[largeX] = asymp_macdo(dprf0[largeX])
    smallX = np.logical_not(largeX)
    res[smallX] = -macdo_x56(dprf0[smallX], k=k)
    return (k1 * L0**(5./3)) * res 





def DPHI(x,y,L0):
    """ 
    dphi = DPHI(x,y,L0) * r0**(-5./3)

   Computes the phase structure function for a separation (x,y).
   The r0 is not taken into account : the final result of DPHI(x,y,L0)
   has to be scaled with r0**-5/3, with r0 expressed in meters, to get
   the right value.
    """

    r = np.sqrt(x**2+y**2)

    """  BEFORE ...... when rod conan did not exist
    fracDim = 5./3.    # Can vary fracDim for those who do not believe in Kolmogorov...
    r53 = r**(fracDim)
    return 6.88*r53
    """ 
    # With L0 ......
    return rodconan(r, L0)


'''
███████╗██╗      █████╗ ██╗   ██╗██╗███╗   ██╗ ██████╗
██╔════╝██║     ██╔══██╗██║   ██║██║████╗  ██║██╔════╝
███████╗██║     ███████║██║   ██║██║██╔██╗ ██║██║  ███╗
╚════██║██║     ██╔══██║╚██╗ ██╔╝██║██║╚██╗██║██║   ██║
███████║███████╗██║  ██║ ╚████╔╝ ██║██║ ╚████║╚██████╔╝
╚══════╝╚══════╝╚═╝  ╚═╝  ╚═══╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝


## ---(Fri May 24 09:33:29 2019)---
runfile('/home/egendron/codes/ADOPT/projects/atlas_ALPAO4K/dm4K.py', wdir='/home/egendron/codes/ADOPT/projects/atlas_ALPAO4K')
nx, ny = 64, 64
Nactu = 3228
x, y, mask = generateActuCoordinates(nx, ny, Nactu)

runfile('/home/egendron/codes/hraa/tools/modesDM.py', wdir='/home/egendron/codes/hraa/tools')
x.shape
x
plt.scatter(x,y,s=1)

selRad=20
xc, yc, xs, ys, idc, ids = selectSlavedList(x,y,selRad)
nc = xc.shape[0]
ns = xs.shape[0]

Ccc = -computeCrossDistance2Matrix(xc, yc, xc, yc)**(5./6)
Csc = -computeCrossDistance2Matrix(xs, ys, xc, yc)**(5./6)
Ccc_1 = np.linalg.pinv(Ccc, rcond=1e-5)
K = Csc.dot(Ccc_1)


Etmp = np.concatenate((np.eye(nc), K), axis=0)
E = np.zeros_like(Etmp)


runfile('/home/egendron/codes/hraa/tools/modesDM.py', wdir='/home/egendron/codes/hraa/tools')
xc, yc, xs, ys, idc, ids = selectSlavedList(x,y,selRad)
E[idc, :] = Etmp[0:nc, :]
E[ids, :] = Etmp[nc:, :]
pl4K(mask, np.sum(E,axis=1))
pl4K(mask, E[:,0])
pl4K(mask, E[:,1])
pl4K(mask, E[:,2])
pl4K(mask, E.dot(np.arange(nc)))
runfile('/home/egendron/codes/hraa/tools/modesDM.py', wdir='/home/egendron/codes/hraa/tools')
l, U = np.linalg.eigh(Ccc)
U.shape
Ufull = np.zeros((3228,1264))
Ufull[idc,:] = U
pl4K(mask, Ufull[:,0])
pl4K(mask, Ufull[:,1])
pl4K(mask, Ufull[:,2])
U = U[:,::-1]
pl4K(mask, Ufull[:,2])
Ufull = Ufull[:,::-1]
pl4K(mask, Ufull[:,2])
pl4K(mask, Ufull[:,3])
pl4K(mask, Ufull[:,4])
pl4K(mask, Ufull[:,5])
Uext = E.dot(U)
pl4K(mask, Uext[:,5])
k=2
pl4K(mask, Uext[:,k])
pl4K(mask, Uext[:,k]); print(k); k+=1
k=55
pl4K(mask, Uext[:,k]); print(k); k+=1
l
Ccc_1 = np.linalg.pinv(Ccc, rcond=1e-5)
K = Csc.dot(Ccc_1)
Etmp = np.concatenate((np.eye(nc), K), axis=0)
E = np.zeros_like(Etmp)
E[idc, :] = Etmp[0:nc, :]
E[ids, :] = Etmp[nc:, :]
Uext = E.dot(U)
k=55
pl4K(mask, Uext[:,k]); print(k); k+=1
k+=10
pl4K(mask, Uext[:,k]); print(k); k+=10
np.diag(Ccc)
np.diag(Ccc).shape
Ccc_1 = np.linalg.inv(Ccc + np.eye(nc))
K = Csc.dot(Ccc_1)
Etmp = np.concatenate((np.eye(nc), K), axis=0)
E = np.zeros_like(Etmp)
E[idc, :] = Etmp[0:nc, :]
E[ids, :] = Etmp[nc:, :]
Uext = E.dot(U)
pl4K(mask, Uext[:,k]); print(k); k+=10
k=1
pl4K(mask, Uext[:,k]); print(k); k+=10
pl4K(mask, Ufull[:,k]); print(k); k+=10
k=0
pl4K(mask, Ufull[:,k]); print(k); k+=1
k=0
pl4K(mask, Uext[:,k]); print(k); k+=1
k+=10
pl4K(mask, Uext[:,k]); print(k); k+=1
pl4K(mask, Uext[:,k]); print(k); k+=10
pl4K(mask, E[:,0])
pl4K(mask, E[:,1])
pl4K(mask, E[:,2])
pl4K(mask, E[:,3])
pl4K(mask, E[:,333])
quit



'''





def computeCrossDistance2Matrix(x1, y1, x2, y2):
    return (x1[:, None] - x2[None, :])**2 + (y1[:, None] - y2[None, :])**2


def rajouteTilt(mat, nt):
    """
    Replaces the first two modes of the basis by pure tip-tilts. Technically,
    this adds 2 rows at the beginning of matrix mat with an identity matrix
    in [0:2, 0:2].

    Parameters
    ----------
    mat : 2D np.array. Modal basis. First 2 modes are supposed to be tip and tilt.
    nt  : int. Number of tilts.

    Returns
    -------
    mat : 2D np.array, modified modal basis.
    """
    # Gestion du tiptilt pour compass
    ni, nj = mat.shape
    mat = np.concatenate((np.zeros((nt, nj)), mat))
    mat[:, :nt] = 0
    mat[:nt, :nt] = np.eye(nt)
    return mat


def edgeModes(x, y, selRadExt=19.0, selRadInt=5.0, x0=None, y0=None,
              filterPiston=True, mergeTilt=True, userSelection=None):
    """
    Compute a modal basis, based only on the actuators that are located in a
    donut of ext radius selRadExt, internal radius selRadInt, and centred
    at (x0, y0).
    The modes are extended/extrapolated on the other actuators using a MMSE,
    kolmogorov-based linear extrapolator.

    Parameters
    ----------
    x, y      : 1D float array. Coordinates of actuators.
    selRadExt : float, optional. External donut radius. The default is 19.0.
    selRadInt : float, optional. Internal donut radius. The default is 5.0.
    x0, y0    : float, optional. Coord of centre of donut. The default is None.
    filterPiston : bool, optional. The default is True.
    mergeTilt : bool, optional. Replaces the first two modes of the basis by
                pure tip-tilts. The default is True.
    userSelection : 1D list of indexes, optional. List of commanded actuators.
    Supersedes all the previous selection parameters. The default is None.

    Returns
    -------
    Ufull : 2D numpy matrix. Modal basis.
    """
    # define where the centre is
    if x0 is None:
        x0 = np.average(x)
    if y0 is None:
        y0 = np.average(x)
    # Re-centre the data
    xx = x - x0
    yy = y - y0
    # distance to centre
    r = np.sqrt(xx * xx + yy * yy)
    # selection mask based on radius
    mskComman = np.logical_and(r < selRadExt, r > selRadInt)
    # list of indexes
    idc = np.where(mskComman)
    Ufull = continuityExtension(x, y, idc)
    return Ufull





def continuityExtension(x, y, userSelection, filterPiston=True, mergeTilt=True):
    """
    Compute a modal basis, based only on the actuators that are selected by
    the index list <userSelection> (called commanded actuators).
    The modes are extended/extrapolated on the other actuators (slave actus)
    using a MMSE, kolmogorov-based linear extrapolator.

    Parameters
    ----------
    x, y          : 1D float array. Coordinates of actuators.
    userSelection : 1D list of indexes. List of commanded actuators.
    filterPiston  : bool, optional. The default is True.
    mergeTilt     : bool, optional. Replaces the first two modes of the basis
                    by pure tip-tilts. The default is True.
    Returns
    -------
    Ufull : 2D numpy matrix. Modal basis.
    """
    # get list of slave actuators
    actu_slaved = np.ones_like(x, dtype=bool)
    idc = userSelection
    actu_slaved[idc] = False
    ids = np.where(actu_slaved)
    # define coords of commanded and slaved actus
    xc = x[idc]
    yc = y[idc]
    xs = x[ids]
    ys = y[ids]
    nc = xc.size  # nombre de commandes
    ns = xs.size  # nbre de slaves
    print(f'Number of commanded actus : {nc}')
    print(f'Number of slaved actus    : {ns}')
    # compute covariance matrix: Ccc covar commnd-command
    Ccc = -computeCrossDistance2Matrix(xc, yc, xc, yc)**(5. / 6)
    # covariance slaved-command
    Csc = -computeCrossDistance2Matrix(xs, ys, xc, yc)**(5. / 6)
    # MMSE for control/extrapolation
    #Ccc_1 = np.linalg.pinv(Ccc, rcond=1e-5)
    Ccc_1 = np.linalg.inv(Ccc + np.eye(nc))
    K = Csc.dot(Ccc_1)
    # Matrix Etmp for prolongation
    Etmp = np.concatenate((np.eye(nc), K), axis=0)
    # mise dans le bon ordre des actionneurs
    E = np.zeros_like(Etmp)
    E[idc, :] = Etmp[0:nc, :]
    E[ids, :] = Etmp[nc:, :]

    # Gestion du piston .. enfin on va essayer ...
    if filterPiston:
        print('Filtering piston out')
        P = computePistonFilteringMatrix(nc)
        Ccc = P.dot(Ccc.dot(P))
        print('Diagonalisation + filtering piston')
        l, U = np.linalg.eigh(Ccc)
        U = U[:, :1:-1]  # Remove piston mode and swap ordering
    else:
        print('Diagonalisation, no filtering of piston')
        l, U = np.linalg.eigh(Ccc)
        U = U[:, ::-1] # swap ordering
        
    # Application of the extension/continuity matrix to the basis U
    Ufull = E.dot(U)
    
    # management of tiptilt
    if mergeTilt==True:
        # Gestion du tiptilt pour compass
        Ufull = rajouteTilt(Ufull, int(2))
    return Ufull


"""
███╗   ███╗ █████╗ ████████╗██████╗ ██╗██╗  ██╗
████╗ ████║██╔══██╗╚══██╔══╝██╔══██╗██║╚██╗██╔╝
██╔████╔██║███████║   ██║   ██████╔╝██║ ╚███╔╝
██║╚██╔╝██║██╔══██║   ██║   ██╔══██╗██║ ██╔██╗
██║ ╚═╝ ██║██║  ██║   ██║   ██║  ██║██║██╔╝ ██╗
╚═╝     ╚═╝╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝╚═╝  ╚═╝

███╗   ███╗ █████╗ ███╗   ██╗██╗██████╗
████╗ ████║██╔══██╗████╗  ██║██║██╔══██╗
██╔████╔██║███████║██╔██╗ ██║██║██████╔╝
██║╚██╔╝██║██╔══██║██║╚██╗██║██║██╔═══╝
██║ ╚═╝ ██║██║  ██║██║ ╚████║██║██║
╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝╚═╝



"""


def insertMatrix(source, dest, lines, cols):
    """
    Insert a matrix within another one.
    <source> is the (little) matrix to be inserted somewhere
    <dest>   is the destination matrix where the source will be inserted
    <lines>  is the list of lines where <source> will be put
    <cols>   is  "    "  " columns  "     "       "   "  "

    Whener required, the destination matrix <dest> will be extended with
    zeros.

    Example
    ttact = [4303, 4304]
    ttmodes = [4301, 4302]
    mat = np.zeros((3,3))
    insertMatrix( np.eye(2), mat, ttact, ttmodes )

    """
    # Transforms lists in numpy arrays
    lines = np.array(lines)
    cols = np.array(cols)
    # searches for limits
    imax = np.max(lines)
    jmax = np.max(cols)
    # compare to size of destination
    idest, jdest = dest.shape
    result = dest.copy()
    if idest <= imax:
        result = np.concatenate((result, np.zeros((imax - idest + 1, jdest))),
                                axis=0)
        idest, jdest = result.shape
    if jdest <= jmax:
        result = np.concatenate((result, np.zeros((idest, jmax - jdest + 1))),
                                axis=1)
        idest, jdest = result.shape
    # copy de la matrice
    for i in range(len(lines)):
        iligne = int(lines[i])
        result[iligne, cols] = source[i, :]
    return result


def insertElem(intList, index, nelem):
    """
    Modifies a list of the indices pointing to elements in an array, when
    a number <nelem> of some new elements are inserted in the array at
    position <index>.

    """
    xList = np.array(intList)
    xList[xList >= index] += nelem
    return list(xList)


def moreLines(mat, nline, follow=None, pos=None):
    """
    Adds a number of <nline> lines to the matrix <mat>.
    These lines can be added at position <pos>.
    By default (pos=None) they are appended at the end of the matrix.

    The argument follow= is a list of lists of indexes of things that
    are in the input matrix. It is modified on output according to the
    matrix extension.

    Example:
        mat = np.ones((10,12))
        iactu = np.arange(10)                # list of actuators
        mat, ipetals, _ = moreLines(mat, 6)  # add 6 petal actuators at the end
        # now add 2 tip and tilt at the beginning
        mat, itt, [iactu, ipetals] = moreLines(mat, 2, pos=0, follow=[iactu, ipetals])

    """
    idest, jdest = mat.shape
    if pos == None:
        pos = idest
    if pos < 0:
        pos = pos + idest

    # defines where are the lines of matrix mat in matrix mat
    xmat = np.arange(idest)
    # defines their new position in the new matrix
    xmat = insertElem(xmat, pos, nline)
    # Pour the coeffs of mat into newmat
    newmat = np.zeros((idest + nline, jdest))
    newmat[xmat, :] = mat

    imore = list(np.arange(nline) + pos)

    if follow != None:
        for i in range(len(follow)):
            follow[i] = insertElem(follow[i], pos, nline)

    return newmat, imore, follow


def moreColumns(mat, ncol, follow=None, pos=None):
    """
    Adds a number of <ncol> columns to the matrix <mat>.
    These columns can be added at position <pos>.
    By default (pos=None) they are appended at the end of the matrix.

    The argument follow= is a list of lists of indexes of things that
    are in the input matrix. It is modified on output according to the
    matrix extension.

    Example:
        mat = np.ones((10,12))
        jactu = np.arange(10)                  # list of actuators
        # add 6 petal actus at the end
        mat, ipetals, _ = moreColumns(mat, 6)
        # now add 2 tip and tilt at the beginning
        mat, itt, [iactu, ipetals] = moreLines(mat, 2, pos=0, follow=[iactu, ipetals])

    """
    idest, jdest = mat.shape
    if pos == None:
        pos = jdest
    if pos < 0:
        pos = pos + jdest

    # defines where are the columns of matrix mat in matrix mat
    xmat = np.arange(jdest)
    # defines their new position in the new matrix
    xmat = insertElem(xmat, pos, ncol)
    # Pour the coeffs of mat into newmat
    newmat = np.zeros((idest, jdest + ncol))
    newmat[:, xmat] = mat

    jmore = list(np.arange(ncol) + pos)

    if follow != None:
        for i in range(len(follow)):
            follow[i] = insertElem(follow[i], pos, ncol)

    return newmat, jmore, follow

def computeMmseReconstructor(dist, ipos_v, ipos_u, L0=None, r0=None, alpha=0, rcond=None):
    """
    delta_pos: matrix of distances between actuactors, meters
    ipos_mmse: actu positions to be mmse-er 
    ipos_keep: actu positions taken as reference 
    L0: outer scale, meters <float>
    r0: Fried parameter, meters <float>
    alpha: factor to lower effect of mmse <float [0,1[>
    """

    if L0 is None:
        L0 = 25    # ou un truc grand
    if r0 is None:
        r0 = 0.144  # un truc médian 

    d_uu = dist[ipos_u, :][:, ipos_u]
    d_uv = dist[ipos_u, :][:, ipos_v]
    infty_cst = rodconan(np.array([1e3]), L0) # converg. value toward infinity
    mult_factor = -0.5 * r0 ** (-5 / 3) * (0.5 / 2 / np.pi) ** 2
    c_uu = mult_factor * (infty_cst + rodconan(d_uu, L0))
    c_uv = mult_factor * (infty_cst + rodconan(d_uv, L0))

    c_uu = c_uu + alpha * np.diag(c_uu) * np.eye(len(c_uu))
    if rcond is not None:
        c_uu_inv = np.linalg.pinv(c_uu, rcond = rcond)
    else:
        c_uu_inv = np.linalg.pinv(c_uu)
    R = np.dot(c_uv.T, c_uu_inv)	# reconstructeur
    return R


def computeMmseMatrix(mat, delta_pos, ipos_mmse, ipos_keep, L0=None, r0=None, alpha=0):
    """
    mat: input matrix to be MMSE-er
    delta_pos: matrix of distances between actuactors, meters
    ipos_mmse: actu positions to be mmse-er 
    ipos_keep: actu positions taken as reference 
    L0: outer scale, meters <float>
    r0: Fried parameter, meters <float>
    alpha: factor to lower effect of mmse <float [0,1[>

    Other comments:
    To compute delta_pos:
    delta_pos = np.sqrt((x[:, None] - x[None, :])**2 + (y[:, None] - y[None, :])**2)    # distances entre actus [mètres]
    if alpha = 0: nominal mmse
    if alpha = 1: mmse effect is null

    Example:
    To extend a basis to the actuators of the ring:
    Br = computeMmseMatrix(B, delta_pos, ipos_ring, ipos_pup, L0=25, r0=ao.turbu.r0, alpha=0.2)
    """

    R = computeMmseReconstructor(delta_pos, ipos_mmse, ipos_keep, L0=L0, r0=r0, alpha=alpha)
    matR= mat.copy()
    matR[ipos_mmse, :] = R.dot(mat[ipos_keep,:])

    return matR

