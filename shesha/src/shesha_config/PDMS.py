#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Created on 12 juil. 2017

@author: vdeo
'''

import numpy as np
from . import shesha_constants as const


#################################################
# P-Class (parametres) Param_dm
#################################################
class Param_dm:

    def __init__(self, debug=False):

        # DM properties
        self.nact = 0  # DM number of actuators
        self.alt = 0.0  # DM conjugation altitude
        self.thresh = 0.0  # Threshold on response for selection
        self.coupling = 0.2  # Actuator coupling (< .3)
        self.hyst = 0.  # DM hysteresis (< 1.)
        self.gain = 1.0  # Actuator gains
        self.pupoffset = np.array([0, 0])
        # Global offset in pupil (x,y) of the whole actuator pattern

        self.unitpervolt = 0.01
        # Influence function sensitivity in unit/volt. Optional [0.01]
        # Stackarray: mic/volt, Tip-tilt: arcsec/volt.
        self.pushforimat = 1.  # nominal voltage for imat

        # Margins for actuator selection
        self.margin = 0.  # LEGACY - outer margin
        self.margin_out = 0.  # outer margin (pitches) from pupil diameter
        self.margin_in = 0.  # inner margin (pitches) from central obstruction
        self.pzt_extent = 5.  # Extent of pzt DM (pitches)

        # KL DM storage structure
        self.nkl = 0  # Number of KL for KL dm
        self._klbas = KL_basis_struct()

        # Hidden variable safe-typed in shesha_constants
        self._type_dm = b''  # Private storage of type_dm
        self._type_pattern = b''  # Private storage of type_pattern
        self._influType = b''  # Private storage of influType
        self._kl_type = b''  # Private storage for KL type

        # HDF5 storage management
        self.file_influ_hdf5 = b''  # Filename for influ hdf5 file
        self.center_name = b''  # Center name in hdf5
        self.cube_name = b''  # Influence function cube name in hdf5
        self.x_name = b''  # x coord name in hdf5
        self.y_name = b''  # y coord name in hdf5
        self.influ_res = b''  # influence resolution name in hdf5
        self.diam_dm = b''  # name for dm diameter
        self.diam_dm_proj = b''  # name for dm diameter in pupil plane

        # PXD cleanup
        # internal kwrd
        self._pitch = None
        """ inter-actuator space in pixels"""
        self._ntotact = None
        """ total number of actuators"""
        self._influsize = None
        """ total number of actuators"""
        self._n1 = None
        """ position of leftmost pixel in largest support"""
        self._n2 = None
        """ position of rightmost pixel in largest support"""
        self._puppixoffset = None
        self._influ = None
        """ influence functions"""
        self._xpos = None
        """ x positions of influ functions"""
        self._ypos = None
        """ y positions of influ functions"""
        self._i1 = None
        self._j1 = None
        self._pupil = None
        """ pupil mask for this dm"""
        self._com = None
        """ current command"""
        self._influpos = None
        self._ninflu = None
        """Influence functions"""
        self._influstart = None  # np.ndarray - Influence function handling
        self._influkernel = None  # np.ndarray convolution kernel for comp_dmshape

    # type_dm : Type of DM
    @property
    def type_dm(self):
        return self._type_dm

    @type_dm.setter
    def type_dm(self, value):
        self._type_dm = const.check_enum(const.DmType, value)

    # type_pattern : type of pattern for pzt DM
    @property
    def type_pattern(self):
        return self._type_pattern

    @type_pattern.setter
    def type_pattern(self, value):
        self._type_pattern = const.check_enum(const.PatternType, value)

    # influType : influence function type for pzt DM
    @property
    def influType(self):
        return self._influType

    @influType.setter
    def influType(self, value):
        self._influType = const.check_enum(const.InfluType, value)

    # kl_type : KL type for KL DM
    @property
    def kl_type(self):
        return self._kl_type

    @kl_type.setter
    def kl_type(self, value):
        self.kl_type = const.check_enum(const.KLType, value)

    def set_pzt_extent(self, p):
        """Set extent of pzt dm in pich unit default = 5
        :param p: (int) : extent pzt dm
        """
        self.pzt_extent = p

    def set_influType(self, t):
        """Set the influence function type for pzt DM
        :param t: (str) : centroider type
        """
        self.influType = t

    def set_gain(self, g):
        """Set the gain to apply to the actuators of the dm

        :param g: (float) : gain
        """
        self.gain = g

    def set_nkl(self, n):
        """Set the number of KL modes used for computation of covmat in case of minimum variance controller

        :param n: (long) : number of KL modes
        """
        self.nkl = n

    def set_kl_type(self, t):
        """Set the type of KL used for computation
        :param t: (string) : KL types : kolmo or karman
        """
        self.kl_type = t

    def set_type(self, t):
        """set the dm type

        :param t: (str) : type of dm
        """
        self.type_dm = t

    def set_pattern(self, t):
        """set the pattern type

        :param t: (str) : type of pattern
        """
        self.type_pattern = t

    def set_file_influ_hdf5(self, f):
        """set the name of hdf5 influence file

        :param filename: (str) : Hdf5 file influence name
        """
        self.file_influ_hdf5 = f

    def set_center_name(self, f):
        """set the name of hdf5 influence file

        :param filename: (str) : Hdf5 file influence name
        """
        self.center_name = f

    def set_cube_name(self, cubename):
        """set the name of influence cube in hdf5

        :param cubename: (str) : name of influence cube
        """
        self.cube_name = cubename

    def set_x_name(self, xname):
        """set the name of x coord of influence fonction in file

        :param t: (str) : name of x coord of influence
        """
        self.x_name = xname

    def set_y_name(self, yname):
        """set the name of y coord of influence fonction in file

        :param yname: (str) : name of y coord of influence
        """
        self.y_name = yname

    def set_influ_res(self, res):
        """set the name of influence fonction resolution in file

        :param res: (str) : name of resoltion (meter/pixel) of influence
        """
        self.influ_res = res

    def set_diam_dm(self, di):
        """set the name of dm diameter in file

        :param di: (str) : name of diameter (meter) dm
        """
        self.diam_dm = di

    def set_diam_dm_proj(self, dp):
        """set the name of dm diameter projet on puille in file

        :param dp: (str) : name of diameter (meter in pupil plan) dm
        """
        self.diam_dm_proj = dp

    def set_nact(self, n):
        """set the number of actuator

        :param n: (long) : number of actuators in the dm
        """
        self.nact = n

    def set_margin(self, n):
        """set the margin for outside actuator select

        :param n: (float) : pupille diametre ratio for actuator select
        """
        self.margin = n

    def set_margin_out(self, n):
        """set the margin for outside actuator select

        :param n: (float) : unit is actuator pitch (+) for extra (-) for intra
        """
        self.margin_out = n

    def set_margin_in(self, n):
        """set the margin for inside actuator select (central obstruction)

        :param n: (float) : unit is actuator pitch (+) for extra (-) for intra
        """
        self.margin_in = n

    def set_alt(self, a):
        """set the conjugaison altitude

        :param a: (float) : conjugaison altitude (im m)
        """
        self.alt = a

    def set_thresh(self, t):
        """set the threshold on response for selection

        :param t: (float) : threshold on response for selection (<1)
        """
        self.thresh = t

    def set_coupling(self, c):
        """set the actuators coupling

        :param c: (float) : actuators coupling (<0.3)
        """
        self.coupling = c

    def set_unitpervolt(self, u):
        """set the Influence function sensitivity

        :param u: (float) : Influence function sensitivity in unit/volt
        """
        self.unitpervolt = u

    def set_push4imat(self, p):
        """set the nominal voltage for imat

        :param p: (float) : nominal voltage for imat
        """
        self.push4imat = p

    def set_ntotact(self, n):
        """set the total number of actuators

        :param n: (long) : total number of actuators
        """
        self._ntotact = n

    def set_xpos(self, xpos):
        """Set the x positions of influ functions (lower left corner)

        :param xpos: (np.ndarray[ndim=1,dtype=np.float32_t]) : x positions of influ functions
        """
        self._xpos = xpos.copy()

    def set_ypos(self, ypos):
        """Set the y positions of influ functions (lower left corner)

        :param ypos: (np.ndarray[ndim=1,dtype=np.float32_t]) : y positions of influ functions
        """
        self._ypos = ypos.copy()

    def set_i1(self, i1):
        """TODO doc

        :param i1: (np.ndarray[ndim=1,dtype=np.int32_t]) :
        """
        self._i1 = i1.copy()

    def set_j1(self, j1):
        """TODO doc

        :param j1: (np.ndarray[ndim=1,dtype=np.int32_t]) :
        """
        self._j1 = j1.copy()

    def set_influ(self, influ):
        """Set the influence function

        :param influ: (np.ndarray[ndim=3,dtype=np.float32_t]) : influence function
        """
        self._influ = influ.copy()


class KL_basis_struct:
    """
        Structure for kl dm parameter
        :param kl: (obj) : kl dm parameter
    """
    nr = 0
    npp = 0
    nfunc = 0
    nord = 0
    radp = None
    evals = None
    npo = None
    ordd = None
    rabas = None
    azbas = None
    ncp = 0
    ncmar = 0
    px = None
    py = None
    cr = None
    cp = None
    pincx = None
    pincy = None
    pincw = None
    ap = None
    kers = None
