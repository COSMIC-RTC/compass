cimport numpy as np

include "sutra.pxd"

cpdef long RASC
cpdef bytes shesha_savepath

#################################################
# P-Class (parametres) Param_loop
#################################################
cdef class Param_loop:
    cdef readonly np.ndarray devices
    """ list of GPU devices """
    cdef readonly long niter
    """ number of iterations"""
    cdef readonly float ittime
    """ iteration time (in sec)"""



#################################################
# P-Class (parametres) Param_tel
#################################################
cdef class Param_tel:
  cdef readonly float diam
  """telescope diameter (in meters)."""
  cdef readonly float cobs
  """central obstruction ratio."""
  cdef readonly str type_ap
  """EELT aperture type: "Nominal", "BP1", "BP3", "BP5" (for back-up plan with 1, 3, or 5 missing annulus)."""
  cdef readonly float t_spiders
  """secondary supports ratio."""
  cdef readonly str spiders_type
  """secondary supports type: "four" or "six"."""
  cdef readonly float pupangle
  """rotation angle of pupil."""
  cdef readonly long nbrmissing
  """number of missing segments for EELT pupil (max is 20)."""
  cdef readonly float referr
  """std of reflectivity errors for EELT segments (fraction)."""
  cdef readonly float std_piston
  """std of piston errors for EELT segments  """
  cdef readonly float std_tt
  """ std of tip-tilt errors for EELT segments """
  cdef readonly np.ndarray vect_seg
  """vector for define segments numbers need"""

#################################################
# P-Class (parametres) Param_geom
#################################################
cdef class Param_geom:
    cdef readonly long  ssize
    """linear size of full image (in pixels)."""
    cdef readonly float zenithangle
    """observations zenith angle (in deg)."""

    # internal keywords
    cdef readonly long  pupdiam
    """linear size of total pupil (in pixels)."""
    cdef readonly float cent
    """central point of the simulation."""
    cdef readonly np.ndarray _ipupil  # total pupil (include full guard band)
    cdef readonly np.ndarray _mpupil  # medium pupil (part of the guard band)
    cdef readonly np.ndarray _spupil  # small pupil (without guard band)
    cdef readonly np.ndarray _phase_ab_M1  # Phase aberration in the pupil (small size)
    cdef readonly np.ndarray _phase_ab_M1_m  # Phase aberration in the pupil (medium size)
    cdef readonly np.ndarray _apodizer  # apodizer (same size as small pupil)
    cdef readonly long  _p1  # min x,y for valid points in mpupil
    cdef readonly long  _p2  # max x,y for valid points in mpupil
    cdef readonly long  _n  # linear size of mpupil
    cdef readonly long  _n1  # min x,y for valid points in ipupil
    cdef readonly long  _n2  # max x,y for valid points in ipupil



#################################################
# P-Class (parametres) Param_wfs
#################################################
cdef class Param_wfs:
    cdef readonly str type_wfs
    """type of wfs : "sh" or "pyr"."""
    cdef readonly long nxsub
    """linear number of subaps."""
    cdef readonly long  npix
    """number of pixels per subap."""
    cdef readonly float pixsize
    """pixel size (in arcsec) for a subap."""
    cdef readonly float Lambda
    """observation wavelength (in um) for a subap."""
    cdef readonly float optthroughput
    """wfs global throughput."""
    cdef readonly float fracsub
    """minimal illumination fraction for valid subaps."""
    cdef readonly long openloop
    """1 if in "open-loop" mode (i.e. does not see dm)."""
    cdef readonly float fssize
    """size of field stop in arcsec."""
    cdef readonly str fstop
    """ Fields of the wfs diaphragm shape : "square" or "none" """

    cdef readonly int atmos_seen
    """1 if the WFS sees the atmosphere layers"""
    cdef readonly np.ndarray dms_seen
    """index of dms seen by the WFS"""
    cdef readonly bool error_budget
    """ If True, enable error budget analysis for the simulation"""


# target kwrd
    cdef readonly float xpos
    """guide star x position on sky (in arcsec)."""
    cdef readonly float ypos
    """guide star x position on sky (in arcsec)."""
    cdef readonly float gsalt
    """altitude of guide star (in m) 0 if ngs."""
    cdef readonly float gsmag
    """magnitude of guide star."""
    cdef readonly float zerop
    """detector zero point expressed in ph/m**2/s in the bandwidth of the WFS"""
    cdef readonly float noise
    """desired noise : < 0 = no noise / 0 = photon only / > 0 photon + ron."""

    cdef readonly float kernel  #
    cdef np.ndarray _ftkernel  # (float*)

# lgs only
    cdef readonly float lgsreturnperwatt
    """return per watt factor (high season : 10 ph/cm2/s/W)."""
    cdef readonly float laserpower
    """laser power in W."""
    cdef readonly float lltx
    """x position (in meters) of llt."""
    cdef readonly float llty
    """y position (in meters) of llt."""
    cdef readonly str proftype
    """type of sodium profile "gauss", "exp", etc ..."""
    cdef readonly float beamsize
    """laser beam fwhm on-sky (in arcsec)."""

# misalignment

    cdef readonly float G
    """Magnifying factor"""
    cdef readonly float thetaML
    """ WFS rotation angle in the pupil"""
    cdef readonly float dx
    """ X axis misalignment in pixels"""
    cdef readonly float dy
    """ Y axis misalignment in pixels"""
# internal kwrd
    cdef readonly long  _pdiam
    """pupil diam for a subap (in pixel)"""
    cdef readonly long  _Nfft
    """ array size for fft for a subap (in pixel)"""
    cdef readonly long  _Ntot
    """ total size of hr image for a subap (in pixel)"""
    cdef readonly long  _nrebin
    """ rebin factor from hr to binned image for a subap"""
    cdef readonly long  _nvalid
    """ number of valid subaps"""

    cdef readonly float   _nphotons
    """ number of photons per subap"""
    cdef readonly float   nphotons4imat
    """ number of photons per subap used for doing imat"""
    cdef readonly float   _subapd
    """ subap diameter (m)"""
    cdef readonly np.ndarray _fluxPerSub
    """ fraction of nphotons per subap"""
    cdef readonly float   _qpixsize
    """ quantum pixel size for the simulation"""

    cdef readonly np.ndarray _istart
    """ (int*) x start indexes for cutting phase screens"""
    cdef readonly np.ndarray _jstart
    """ (int*) y start indexes for cutting phase screens"""
    # cdef np.ndarray _validsubs    # (i,j) indices of valid subaps
    cdef readonly np.ndarray _validsubsx
    """(int*) indices of valid subaps along axis x"""
    cdef readonly np.ndarray _validsubsy
    """(int*) indices of valid subaps along axis y"""
    cdef readonly np.ndarray _isvalid
    """(int*) array of 0/1 for valid subaps"""
    cdef readonly np.ndarray _phasemap
    """(int*) array of pixels transform from phase screen into subaps phase screens"""
    cdef readonly np.ndarray _hrmap
    """(int*) array of pixels transform from minimal FoV image to (in case type is sh or geo)"""
    cdef readonly np.ndarray _sincar
    """(float*) array of pixels transform from minimal FoV image to (in case type is "pyr" or "roof")"""
                         # full FoV image (array of 0 if the same)
    cdef readonly np.ndarray _binmap
    """(int*) array of pixels transform from full FoV hr images to binned images"""
    cdef readonly np.ndarray _halfxy
    """(float*) phase offset for 1/2 pixel shift in (x,y)"""

    cdef readonly np.ndarray _submask
    """(float*) fieldstop for each subap"""

    cdef readonly np.ndarray _lgskern
    """ lgs kernels for each subap"""
    cdef readonly np.ndarray _profna
    """ sodium profile"""
    cdef readonly np.ndarray _altna
    """ corresponding altitude"""
    cdef readonly np.ndarray _prof1d
    """ hr profile"""
    cdef readonly np.ndarray _profcum
    """ hr profile cumulated"""
    cdef readonly np.ndarray _beam
    """ 1d beam function"""
    cdef readonly np.ndarray _ftbeam
    """ 1d beam function fft"""
    cdef readonly np.ndarray _azimuth
    """ angles of rotation for each spot"""

# pyramid-nly kwrds
    cdef readonly float pyr_ampl
    """pyramid wfs modulation amplitude radius [arcsec]."""
    cdef readonly long pyr_npts
    """total number of point along modulation circle [unitless]."""
    cdef np.ndarray pyr_pos
    """positions for modulation, overwrites ampl and npts [arcsec]"""
    cdef readonly str pyr_loc
    """Location of modulation, before/after the field stop.
    valid value are "before" or "after" (default "after")."""
    cdef readonly str pyrtype
    """Type of pyramid, either 0 for "Pyramid" or 1 for "RoofPrism"."""
    cdef readonly long pyr_pup_sep
    """Pyramid pupil separation. (default: long(wfs.nxsub))"""

# pyramid internal kwrds
    cdef readonly np.ndarray _pyr_offsets  # (float*)
    cdef readonly np.ndarray _pyr_cx  # (float*)
    cdef readonly np.ndarray _pyr_cy  # (float*)


    # cdef make_lgs_prof1d(self, Param_tel p_tel,
    #        np.ndarray[dtype=np.float32_t] prof,
    #        np.ndarray[dtype=np.float32_t] h,
    #        float beam, bytes center=?)



#################################################
# P-Class (parametres) Param_atmos
#################################################
cdef class Param_atmos:

    cdef readonly long nscreens
    """number of turbulent layers."""
    cdef readonly float r0
    """global r0."""
    cdef readonly float pupixsize
    """pupil pixel size (in meters)."""
    cdef readonly np.ndarray L0
    """L0 per layers in meters."""
    cdef readonly np.ndarray dim_screens
    """linear size of phase screens."""
    cdef readonly np.ndarray alt
    """altitudes of each layer."""
    cdef readonly np.ndarray winddir
    """wind directions of each layer."""
    cdef readonly np.ndarray windspeed
    """wind speeds of each layer."""
    cdef readonly np.ndarray frac
    """fraction of r0 for each layer."""
    cdef readonly np.ndarray deltax
    """x translation speed (in pix / iteration) for each layer."""
    cdef readonly np.ndarray deltay
    """y translation speed (in pix / iteration) for each layer."""
    cdef readonly np.ndarray seeds


#

cdef class Klbas:
    """ structure for kl dm parameter
    :param kl: (obj) : kl dm parameter
    """

    cdef readonly long nr
    cdef readonly long npp
    cdef readonly long nfunc
    cdef readonly long nord
    cdef readonly np.ndarray radp
    cdef readonly np.ndarray evals
    cdef readonly np.ndarray npo
    cdef readonly np.ndarray ordd
    cdef readonly np.ndarray rabas
    cdef readonly np.ndarray azbas
    cdef readonly long ncp
    cdef readonly long ncmar
    cdef readonly np.ndarray px
    cdef readonly np.ndarray py
    cdef readonly np.ndarray cr
    cdef readonly np.ndarray cp
    cdef readonly np.ndarray pincx
    cdef readonly np.ndarray pincy
    cdef readonly np.ndarray pincw
    cdef readonly np.ndarray ap
    cdef readonly np.ndarray kers

#################################################
# P-Class (parametres) Param_dm
#################################################
cdef class Param_dm:
  cdef readonly bytes  influType
  """ type of influence function for pzt dms"""

  cdef readonly bytes  type_dm
  """ type of dm"""
  cdef readonly bytes  type_pattern
  """ type of pattern"""
  cdef readonly bytes  file_influ_hdf5
  """ filename for influ hdf5 file"""
  cdef readonly bytes  center_name
  """ filename for influ hdf5 file"""
  cdef readonly bytes  cube_name
  """ name for influence cube in hdf5"""
  cdef readonly bytes  x_name
  """ name for x coord of influence"""
  cdef readonly bytes  y_name
  """ name for y coord of influence"""
  cdef readonly bytes influ_res
  """ name for influence resolution"""
  cdef readonly bytes diam_dm
  """ name for diameter dm"""
  cdef readonly bytes diam_dm_proj
  """ name for diameter dm in pupil plan"""
  cdef readonly long    nact
  """ number of actuators in dm """
  cdef readonly float   alt
  """ conjugaison altitude (im m)"""
  cdef readonly float   thresh
  """ threshold on response for selection (<1)"""
  cdef readonly float   coupling
  """ actuators coupling (<0.3)"""
  cdef readonly float   hyst
  """ actuators hysteresis (<1.)"""
  cdef readonly float   margin
  """outside margin for actuator select"""
  cdef readonly float   margin_out
  """outside margin for actuator select"""

  cdef readonly float   margin_in
  """inside margin for actuator select"""
  cdef readonly float   gain
  """gain to apply to the actuator of the dm"""

  cdef readonly np.ndarray   pupoffset
  """(2)"""
  """ global offset in pupil of whole actuator pattern [m]"""
  cdef readonly float   unitpervolt
  """ Influence function sensitivity in unit/volt. Optional [0.01]
      Stackarray: mic/volt, Tip-tilt: arcsec/volt."""
  cdef readonly float   push4imat
  """ nominal voltage for imat"""
  cdef readonly long    nkl
  """ number of kl modes"""
  cdef readonly bytes   kl_type
  """ kl_type : kolmo or karman"""

  # internal kwrd
  cdef readonly long    _pitch
  """ inter-actuator space in pixels"""
  cdef readonly long    _ntotact
  """ total number of actuators"""
  cdef readonly long    _influsize
  """ total number of actuators"""
  cdef readonly long    _n1
  """ position of leftmost pixel in largest support"""
  cdef readonly long    _n2
  """ position of rightmost pixel in largest support"""
  cdef readonly long    _puppixoffset[2]
  cdef readonly np.ndarray _influ
  """ influence functions"""
  cdef readonly np.ndarray _xpos
  """ x positions of influ functions"""
  cdef readonly np.ndarray _ypos
  """ y positions of influ functions"""
  cdef readonly np.ndarray _i1
  cdef readonly np.ndarray _j1
  cdef readonly np.ndarray _pupil
  """ pupil mask for this dm"""
  cdef readonly np.ndarray _com
  """ current command"""
  cdef readonly np.ndarray _influpos
  cdef readonly np.ndarray _ninflu
  """Influence functions"""
  cdef readonly np.ndarray _influstart
  cdef readonly np.ndarray _influkernel  # Array to convolution kernel for comp_dmshape

  cdef readonly Klbas _klbas
  """ np.ndarray to a kl struct"""


  # cdef set_xpos(self,np.ndarray[ndim=1,dtype=np.float32_t] xpos)
  # cdef set_ypos(self,np.ndarray[ndim=1,dtype=np.float32_t] ypos)



#################################################
# P-Class (parametres) Param_target
#################################################
cdef class Param_target:
    cdef readonly int ntargets
    """number of targets"""
    cdef readonly int apod
    """boolean for apodizer"""
    cdef readonly np.ndarray Lambda
    """observation wavelength for each target"""
    cdef readonly np.ndarray xpos
    """x positions on sky (in arcsec) for each target"""
    cdef readonly np.ndarray ypos
    """y positions on sky (in arcsec) for each target"""
    cdef readonly np.ndarray mag
    """magnitude for each target"""
    cdef readonly float zerop
    """target flux for magnitude 0"""
    cdef readonly np.ndarray  dms_seen
    """index of dms seen by the target"""



#################################################
# P-Class (parametres) Param_rtc
#################################################
cdef class Param_rtc:
    cdef readonly long    nwfs
    # number of wfs
    cdef readonly list centroiders
    # an array of centroiders
    cdef readonly list controllers
    # an array of controllers



#################################################
# P-Class (parametres) Param_centroider
#################################################
cdef class Param_centroider:
    cdef readonly long    nwfs
    """ index of wfs in y_wfs structure on which we want to do centroiding"""
    cdef readonly bytes  type_centro
    """ type of centroiding cog, tcog, bpcog, wcog, corr"""
    cdef readonly bytes  type_fct
    """ type of ref function gauss, file, model"""
    cdef readonly np.ndarray weights
    """ optional reference function(s) used for centroiding"""
    cdef readonly long    nmax
    """ number of brightest pixels"""
    cdef readonly float   thresh
    """Threshold"""
    cdef readonly float   width
    """ width of the Gaussian"""
    cdef readonly long    sizex
    """ x-size for inter mat (corr)"""
    cdef readonly long    sizey
    """ x-size for inter mat (corr)"""
    cdef readonly np.ndarray interpmat
    """ optional reference function(s) used for corr centroiding"""
    cdef readonly int method
    """ optional method used in the pyrhr centroider (0:local, 1:global)"""


#################################################
# P-Class (parametres) Param_controller
#################################################
cdef class Param_controller:
    cdef readonly bytes  type_control
    """ type of controller"""
    cdef readonly int kl_imat
    """ set imat kl on-off"""
    cdef readonly np.ndarray klgain
    """ gain for kl mod in imat kl """
    cdef readonly np.ndarray nwfs
    """ index of wfss in controller"""
    cdef readonly np.ndarray nvalid
    """ number of valid subaps per wfs"""
    cdef readonly np.ndarray ndm
    """ index of dms in controller"""
    cdef readonly np.ndarray nactu
    """ number of controled actuator per dm"""
    cdef readonly np.ndarray imat
    """ full interaction matrix"""
    cdef readonly np.ndarray cmat
    """ full control matrix"""
    cdef readonly float   maxcond
    """ max condition number"""
    cdef readonly float   TTcond
    """ tiptilt condition number for cmat filtering with mv controller"""
    cdef readonly float    delay
    """ loop delay [frames]"""
    cdef readonly float   gain
    """ loop gain """
    cdef readonly long    nkl
    cdef readonly long    cured_ndivs
    """ subdivision levels in cured"""
    cdef readonly int     modopti
    """ Flag for modal optimization"""
    cdef readonly int     nrec
    """ Number of sample of open loop slopes for modal optimization computation"""
    cdef readonly int     nmodes
    """ Number of modes for M2V matrix (modal optimization)"""
    cdef readonly float     gmin
    """ Minimum gain for modal optimization"""
    cdef readonly float     gmax
    """ Maximum gain for modal optimization"""
    cdef readonly int     ngain
    """ Number of tested gains"""


cpdef  indices(int dim1, int dim2=?)
cpdef rotate3d(np.ndarray[ndim=3, dtype=np.float32_t] im,
              np.ndarray[ndim=1, dtype=np.float32_t] ang,
              float cx=?, float cy=?, float zoom=?)
cpdef rotate(np.ndarray[ndim=3, dtype=np.float32_t] im,
            float ang, float cx=?, float cy=?, float zoom=?)
cpdef makegaussian(int size, float fwhm, int xc=?, int yc=?, int norm=?)
