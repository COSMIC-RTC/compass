from chakra cimport *

cimport numpy as np

from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from libc.stdint cimport uintptr_t


cdef np.float32_t dtor = (np.pi/180.)


ctypedef enum cudaChannelFormatDesc:
    f
    w
    x
    y
    z


#################################################
# C-Class sutra_phase
#################################################
cdef extern from "sutra_phase.h":
    cdef cppclass sutra_phase:
        carma_obj[float] *d_screen
        long screen_size
        float *zerCoeff
        carma_obj[float] *zernikes
        carma_obj[float] *mat


cdef extern from "sutra_turbu.h":
#################################################
# C-Class sutra_t_screen
#################################################
    cdef cppclass sutra_tscreen:
        int device # The device # 
        sutra_phase *d_tscreen # The phase screen   
        long screen_size # size of phase screens
        float amplitude # amplitude for extrusion (r0^-5/6)
        float altitude
        float windspeed
        float winddir
        float deltax # number of rows to extrude per iteration
        float deltay # number of lines to extrude per iteration
        #internal
        float accumx
        float accumy

        float norm_vk
        bool vk_on
        carma_context *current_context

        carma_obj[float] *d_A # A matrix for extrusion
        carma_obj[float] *d_B # B matrix for extrusion
        carma_obj[unsigned int] *d_istencilx # stencil for column extrusion
        carma_obj[unsigned int] *d_istencily # stencil for row extrusion
        carma_obj[float] *d_z # tmp array for extrusion process
        carma_obj[float] *d_noise # tmp array containing random numbers
        carma_obj[float] *d_ytmp

        sutra_tscreen(carma_context *context, long size, long size2, float amplitude,
                      float altitude, float windspeed, float winddir, float deltax,
                      float deltay, int device)
        #sutra_tscreen(const sutra_tscreen& tscreen)

        int init_screen(float *h_A, float *h_B, unsigned int *h_istencilx,
                        unsigned int *h_istencily, int seed)
        int extrude(int dir)
        int init_vk(int seed, int pupd)
        int generate_vk(float l0, int nalias)


#################################################
# C-Class sutra_atmos
#################################################
    cdef cppclass sutra_atmos:
        int nscreens
        map.map[float, sutra_tscreen *] d_screens
        float r0
        carma_obj[float] *d_pupil
        carma_context *current_context

        sutra_atmos(carma_context *context, int nscreens, float *r0, long *size,
                    long *size2, float *altitude, float *windspeed, float *winddir,
                    float *deltax, float *deltay, float *pupil, int device)
        init_screen(float alt, float *h_A, float *h_B, unsigned int *h_istencilx,
                    unsigned int *h_istencily, int seed)
        int move_atmos()



cdef extern from "sutra_target.h":
#################################################
# C-Class sutra_source
#################################################
    cdef cppclass sutra_source:
        
        int device #   device # 
        float tposx #   x position of target on the sky  
        float tposy #   y position of target on the sky
        long npos #   number of points in the pupil
        float mag #   brightness of target
        float Lambda "lambda" #   imaging lambda
        float zp #   imaging zero point
        float scale #   phase scale
        bool lgs #   flag for lgs
        char* type #   type of source : target / wfs
        int block_size #   optimum block size of device

        float strehl_se #   short exposure strehl
        float strehl_le #   long exposure strehl
        float ref_strehl #   reference for strehl computation
        int strehl_counter #   counter for le strehl computation
        float phase_var #   current phase variance in the pupil
        float phase_var_avg #   average phase variance in the pupil
        int phase_var_count #   counter for average phase variance in the pupil

        sutra_phase *d_phase #   phase for this target
        # INTRO PHASE INSTRU
        #sutra_phase *d_phase_instru;
        #
        carma_host_obj[float] *phase_telemetry #   
        # LGS unkwon for now   sutra_lgs *d_lgs #   the lgs object
        carma_obj[float] *object #   the object intensity map
        carma_obj[float] *d_pupil #   the pupil mask
        carma_obj[float] *d_image #   the resulting image for this target
        carma_obj[float] *d_leimage #   the long exposure image for this target
        carma_obj[cuFloatComplex] *d_amplipup #   the complex amplitude in the pupil plane
        carma_obj[float] *d_phasepts #   the valid phase points in the pupil (target only)
        carma_obj[int] *d_wherephase #   positions of valid phase points in the pupil (target only)
        #type_screen unknown cf sutra_DM
        # map[type_screen, float] xoff #   x reference for raytracing
        # map[type_screen, float] yoff #   y reference for raytracing
        carma_context *current_context
        
        int add_layer(char* l_type, float alt, float xoff, float yoff)
        int init_strehlmeter()
        int raytrace(sutra_atmos *atmos)
        int comp_image(int puponly)
        int comp_strehl()
        
    cdef cppclass sutra_target:
        int ntargets
        vector.vector[sutra_source *] d_targets

        sutra_target(carma_context *context, int ntargets, float *xpos, float *ypos,
            float *Lambda, float *mag, long *sizes, float *pupil, int device)
        # not implemented
        #int get_phase(int ntarget, float *dest)


    cdef int target_texraytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
            int Ny, float xoff, float yoff, int Ntot, cudaChannelFormatDesc channelDesc,
            carma_device *device)

    cdef int target_raytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
            float xoff, float yoff, int block_size)

    cdef int target_lgs_raytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
            float xoff, float yoff, float delta, int block_size)

    cdef int target_raytrace_async(carma_streams streams, float *d_odata, float *d_idata,
            int nx, int ny, int Nx, float xoff, float yoff, int block_size)

    cdef int target_raytrace_async(carma_host_obj[float] *phase_telemetry, float *d_odata,
            float *d_idata, int nx, int ny, int Nx, float xoff, float yoff, int block_size)

    cdef int fft_goodsize(long size)



#################################################
# P-Class atmos
#################################################
cdef class atmos:
    cdef sutra_atmos *s_a
    cdef carma_context context


#################################################
# P-Class target
#################################################
cdef class target:
    cdef sutra_target *target
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
    cdef int device
    cdef chakra_context context


#################################################
# P-Class (parametres) param_atmos
#################################################
cdef class param_atmos:

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


#################################################
# P-Class (parametres) param_geom
#################################################
cdef class param_geom:
    cdef readonly long  ssize
    """linear size of full image (in pixels)."""
    cdef readonly float zenithangle
    """observations zenith angle (in deg)."""
  
    # internal keywords
    cdef readonly long  pupdiam
    """linear size of total pupil (in pixels)."""
    cdef readonly float cent
    """central point of the simulation."""
    cdef np.ndarray _ipupil   # total pupil (include full guard band)
    cdef np.ndarray _mpupil   # medium pupil (part of the guard band)
    cdef np.ndarray _spupil   # small pupil (without guard band)
    cdef np.ndarray _apodizer # apodizer (same size as small pupil)
    cdef long  _p1         # min x,y for valid points in mpupil
    cdef long  _p2         # max x,y for valid points in mpupil
    cdef long  _n          # linear size of mpupil
    cdef long  _n1         # max x,y for valid points in ipupil
    cdef long  _n2         # min x,y for valid points in ipupil


#################################################
# P-Class (parametres) param_tel
#################################################
cdef class param_tel:
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


#################################################
# P-Class (parametres) param_target
#################################################
cdef class param_target:
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

#################################################
# P-Class (parametres) param_wfs
#################################################
cdef class param_wfs:
    cdef public str typ
    """type of wfs : "sh" or "pyr"."""
    cdef public long nxsub
    """linear number of subaps."""
    cdef public long  npix
    """number of pixels per subap."""
    cdef public float pixsize
    """pixel size (in arcsec) for a subap."""
    cdef public float Lambda
    """observation wavelength (in µm) for a subap."""
    cdef public float optthroughput
    """wfs global throughput."""
    cdef public float fracsub
    """minimal illumination fraction for valid subaps."""
    cdef public long openloop
    """1 if in "open-loop" mode (i.e. does not see dm)."""
    cdef public float fssize
    """size of field stop in arcsec."""
    cdef public str fstop
    """size of field stop in arcsec."""

    #TODO 
    #int atmos_seen;        // 1 if the WFS sees the atmosphere layers
    #pointer dms_seen;      // index of dms seen by the WFS


#target kwrd
    cdef public float xpos
    """guide star x position on sky (in arcsec)."""
    cdef public float ypos
    """guide star x position on sky (in arcsec)."""
    cdef public float gsalt
    """altitude of guide star (in m) 0 if ngs."""
    cdef public float gsmag
    """magnitude of guide star."""
    cdef public float zerop
    """detector zero point."""
    cdef public float noise
    """desired noise : < 0 = no noise / 0 = photon only / > 0 photon + ron."""

    cdef public float kernel      # 
    cdef float* _ftkernel #

# lgs only
    cdef public float lgsreturnperwatt
    """return per watt factor (high season : 10 ph/cm2/s/W)."""
    cdef public float laserpower
    """laser power in W."""
    cdef public float lltx
    """x position (in meters) of llt."""
    cdef public float llty
    """y position (in meters) of llt."""
    cdef public str proftype
    """type of sodium profile "gauss", "exp", etc ..."""
    cdef public float beamsize
    """laser beam fwhm on-sky (in arcsec)."""

#internal kwrd
    cdef long  _pdiam          # pupil diam for a subap (in pixel)
    cdef long  _Nfft           # array size for fft for a subap (in pixel)
    cdef long  _Ntot           # total size of hr image for a subap (in pixel)
    cdef long  _nrebin         # rebin factor from hr to binned image for a subap 
    cdef long  _nvalid         # number of valid subaps

    cdef float   _nphotons     # number of photons per subap
    cdef float   _subapd       # subap diameter (m)
    cdef float* _fluxPerSub   # fraction of nphotons per subap
    cdef float   _qpixsize     # quantum pixel size for the simulation

    cdef float* _istart       # x start indexes for cutting phase screens 
    cdef float* _jstart       # y start indexes for cutting phase screens 
    cdef float* _validsubs    # (i,j) indices of valid subaps
    cdef float* _isvalid      # array of 0/1 for valid subaps
    cdef float* _phasemap     # array of pixels transform from phase screen into
                         # subaps phase screens
    cdef float* _hrmap        # array of pixels transform from minimal FoV image to
                         # full FoV image (array of 0 if the same)
    cdef float* _binmap       # array of pixels transform from full FoV hr images to
                         # binned images
    cdef float* _halfxy       # phase offset for 1/2 pixel shift in (x,y)

    cdef float* _submask       # fieldstop for each subap

    cdef float* _lgskern      # lgs kernels for each subap
    cdef float* _profna       # sodium profile
    cdef float* _altna        # corresponding altitude
    cdef float* _prof1d       # hr profile
    cdef float* _profcum      # hr profile cumulated
    cdef float* _beam         # 1d beam function
    cdef float* _ftbeam       # 1d beam function fft
    cdef float* _azimuth      # angles of rotation for each spot

# pyramid-nly kwrds
    cdef public float   pyr_ampl
    """pyramid wfs modulation amplitude radius [arcsec]."""
    cdef public long    pyr_npts
    """total number of point along modulation circle [unitless]."""
    cdef float* pyr_pos    # positions for modulation, overwrites ampl and npts [arcsec]
    cdef public str  pyr_loc
    """Location of modulation, before/after the field stop.
    valid value are "before" or "after" (default "after")."""
    cdef public str  pyrtype
    """Type of pyramid, either 0 for "Pyramid" or 1 for "RoofPrism"."""

# pyramid internal kwrds
    cdef float* _pyr_offsets   #
    cdef float* _pyr_cx   #
    cdef float* _pyr_cy   #



#################################################
# P-Class (parametres) param_dm
#################################################
cdef class param_dm:
  cdef str  type
  """ type of dm"""
  cdef readonly long    nact
  """ number of actuators in the diameter """
  cdef readonly float   alt
  """ conjugaison altitude (im m)"""
  cdef readonly float   thresh
  """ threshold on response for selection (<1)"""
  cdef readonly float   coupling
  """ actuators coupling (<0.3)"""
  cdef readonly float   hyst
  """ actuators hysteresis (<1.)"""
  cdef readonly float   margin
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
  
  #internal kwrd
  cdef long    _pitch
  """ inter-actuator space in pixels"""
  cdef long    _ntotact
  """ total number of actuators"""
  cdef long    _influsize
  """ total number of actuators"""
  cdef long    _n1
  """ position of leftmost pixel in largest support"""
  cdef long    _n2
  """ position of rightmost pixel in largest support"""
  cdef long    _puppixoffset[2]
  cdef np.ndarray _influ
  """ influence functions"""
  cdef np.ndarray _xpos
  """ x positions of influ functions"""
  cdef np.ndarray _ypos
  """ y positions of influ functions"""
  cdef np.ndarray _i1
  cdef np.ndarray _j1
  cdef np.ndarray _pupil
  """ pupil mask for this dm"""
  cdef np.ndarray _com
  """ current command"""
  cdef np.ndarray _influpos
  cdef np.ndarray _ninflu 
  cdef np.ndarray _influstart
  cdef np.ndarray _klbas
  """ np.ndarray to a kl struct"""

