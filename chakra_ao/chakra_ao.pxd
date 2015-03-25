from chakra cimport *

import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc, free
from libcpp.cast cimport dynamic_cast
from libcpp.string cimport string
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

    #cdef int fft_goodsize(long size)




cdef extern from "sutra_wfs.h":
#################################################
# C-Class sutra_sensor
#################################################
    cdef cppclass sutra_sensors:
        int device;
        carma_context *current_context;
        size_t nsensors() 
        
        vector[sutra_wfs *] d_wfs;
        map[vector[int],cufftHandle*] campli_plans;
        map[vector[int],cufftHandle*] fttotim_plans;
        map[vector[int],cufftHandle*] ftlgskern_plans;

        carma_obj[cuFloatComplex] *d_camplipup;
        carma_obj[cuFloatComplex] *d_camplifoc;
        carma_obj[cuFloatComplex] *d_fttotim;
        carma_obj[cuFloatComplex] *d_ftlgskern;
        carma_obj[float] *d_lgskern;

        sutra_sensors(carma_context *context, const char** type, int nwfs, long *nxsub,
          long *nvalid, long *npix, long *nphase, long *nrebin, long *nfft,
          long *ntot, long *npup, float *pdiam, float *nphot, int *lgs, int device);
        sutra_sensors(carma_context *context, int nwfs, long *nxsub, long *nvalid,
          long *nphase, long npup, float *pdiam, int device);

        int sensors_initgs(float *xpos, float *ypos, float *Lambda, float *mag,
          long *size, float *noise, long *seed);
        int sensors_initgs(float *xpos, float *ypos, float *Lambda, float *mag,
          long *size, float *noise);
        int sensors_initgs(float *xpos, float *ypos, float *Lambda, float *mag,
          long *size)
        int allocate_buffers()

#################################################
# C-Class sutra_wfs
#################################################
    cdef cppclass sutra_wfs:

        int device
        string type
        long nxsub
        long nvalid
        long npix
        long nrebin
        long nfft
        long ntot
        long npup
        long nphase
        long nmaxhr
        long nffthr
        float subapd
        float nphot
        float noise
        bool lgs
        bool kernconv

        cufftHandle *campli_plan
        cufftHandle *fttotim_plan
        carma_obj[cuFloatComplex] *d_ftkernel
        carma_obj[cuFloatComplex] *d_camplipup
        carma_obj[cuFloatComplex] *d_camplifoc
        carma_obj[cuFloatComplex] *d_fttotim

        carma_obj[float] *d_pupil
        carma_obj[float] *d_bincube
        carma_obj[float] *d_binimg
        carma_obj[float] *d_subsum
        carma_obj[float] *d_offsets
        carma_obj[float] *d_fluxPerSub
        carma_obj[float] *d_sincar
        carma_obj[int] *d_hrmap

        carma_obj[int] *d_isvalid # nxsub x nxsub
        carma_obj[float] *d_slopes

        carma_host_obj[float] *image_telemetry

      # sh only
        carma_obj[int] *d_phasemap
        carma_obj[int] *d_binmap
        carma_obj[int] *d_validsubsx # nvalid
        carma_obj[int] *d_validsubsy # nvalid
        carma_obj[int] *d_istart # nxsub 
        carma_obj[int] *d_jstart # nxsub

      # pyramid only
        carma_obj[float] *d_hrimg
        carma_obj[float] *d_submask
        carma_obj[float] *d_psum
        carma_obj[cuFloatComplex] *d_phalfxy
        carma_obj[cuFloatComplex] *d_poffsets

        carma_host_obj[int] *pyr_cx
        carma_host_obj[int] *pyr_cy

        sutra_source *d_gs

        carma_streams *streams
        int nstreams

        carma_context *current_context

        sutra_wfs(carma_context *context,sutra_sensors *sensors,  const char* type, long nxsub, long nvalid,
          long npix, long nphase, long nrebin, long nfft, long ntot, long npup,
          float pdiam, float nphotons, int lgs, int device)
        sutra_wfs(carma_context *context, long nxsub, long nvalid, long nphase,
          long npup, float pdiam, int device)
        sutra_wfs(const sutra_wfs& wfs)

        int wfs_initarrays(int *phasemap, int *hrmap, int *binmap, float *offsets,
          float *pupil, float *fluxPerSub, int *isvalid, int *validsubsx,
          int *validsubsy, int *istart, int *jstart, cuFloatComplex *kernel)
        int wfs_initarrays(int *phasemap, float *offsets, float *pupil, float *fluxPerSub,
          int *isvalid, int *validsubsx, int *validsubsy, int *istart, int *jstart)
        int wfs_initarrays(cuFloatComplex *halfxy, cuFloatComplex *offsets,
          float *focmask, float *pupil, int *isvalid, int *cx, int *cy,
          float *sincar, int *phasemap, int *validsubsx, int *validsubsy)
        int wfs_initgs(sutra_sensors *sensors, float xpos, float ypos, float Lambda, float mag, long size,
          float noise, long seed)
        int load_kernels(float *lgskern)
        int sensor_trace(sutra_atmos *yatmos)
        # TODO int sensor_trace(sutra_dms *ydm, int rst)
        # TODO int sensor_trace(sutra_atmos *atmos, sutra_dms *ydms)
        int comp_image_tele()
        int fill_binimage()
        int comp_image()
        int slopes_geom(int type, float *slopes)
        int slopes_geom(int type)

        int comp_sh_generic()
        int comp_pyr_generic()
        int comp_roof_generic()


#################################################
# C-Class sutra_wfs_sh
#################################################
cdef extern from "sutra_wfs_sh.h":
    cdef cppclass sutra_wfs_sh(sutra_wfs):
         int  wfs_initarrays(int *phasemap, int *hrmap, int *binmap, float *offsets,
      float *pupil, float *fluxPerSub, int *isvalid, int *validsubsx,
      int *validsubsy, int *istart, int *jstart, cuFloatComplex *kernel)


#################################################
# C-Class sutra_wfs_pyr_roof
#################################################
cdef extern from "sutra_wfs_pyr.h":
    cdef cppclass sutra_wfs_pyr(sutra_wfs):
        int wfs_initarrays(cuFloatComplex *halfxy, cuFloatComplex *offsets,
      float *focmask, float *pupil, int *isvalid, int *cx, int *cy,
      float *sincar, int *phasemap, int *validsubsx, int *validsubsy)

#################################################
# C-Class sutra_wfs_pyr_roof
#################################################
cdef extern from "sutra_wfs_pyr_roof.h":
    cdef cppclass sutra_wfs_pyr_roof(sutra_wfs_pyr):
        sutra_wfs_pyr_roof(const sutra_wfs_pyr_roof& wfs)

#################################################
# C-Class sutra_wfs_pyr_pyr4
#################################################
cdef extern from "sutra_wfs_pyr_pyr4.h":
    cdef cppclass sutra_wfs_pyr_pyr4(sutra_wfs_pyr):
        sutra_wfs_pyr_pyr4(const sutra_wfs_pyr_pyr4& wfs)

#################################################
# C-Class sutra_wfs_geom
#################################################
cdef extern from "sutra_wfs_geom.h":
    cdef cppclass sutra_wfs_geom(sutra_wfs):
        #sutra_wfs_geom(carma_context *context, long nxsub, long nvalid, long nphase,
        #    long npup, float pdiam, int device)
        sutra_wfs_geom(const sutra_wfs_geom& wfs)
        int wfs_initarrays(int *phasemap, float *offsets, float *pupil,
            float *fluxPerSub, int *isvalid, int *validsubsx, int *validsubsy)


cdef extern from *:
    sutra_wfs_geom* dynamic_cast_wfs_geom_ptr "dynamic_cast<sutra_wfs_geom*>" (sutra_wfs*) except NULL
    sutra_wfs_sh* dynamic_cast_wfs_sh_ptr "dynamic_cast<sutra_wfs_sh*>" (sutra_wfs*) except NULL
    sutra_wfs_pyr_pyr4* dynamic_cast_wfs_pyr_pyr4_ptr "dynamic_cast<sutra_wfs_pyr_pyr4*>" (sutra_wfs*) except NULL
    sutra_wfs_pyr_roof* dynamic_cast_wfs_pyr_roof_ptr "dynamic_cast<sutra_wfs_pyr_roof*>" (sutra_wfs*) except NULL

#################################################
# P-Class Atmos
#################################################
cdef class Atmos:
    cdef sutra_atmos *s_a
    cdef chakra_context context
    cdef realinit(self,chakra_context ctxt,int nscreens,
                np.ndarray[dtype=np.float32_t] r0,
                np.ndarray[dtype=np.int64_t] size,
                np.ndarray[dtype=np.float32_t] altitude,
                np.ndarray[dtype=np.float32_t] windspeed,
                np.ndarray[dtype=np.float32_t] winddir,
                np.ndarray[dtype=np.float32_t] deltax,
                np.ndarray[dtype=np.float32_t] deltay,
                np.ndarray[ndim=2,dtype=np.float32_t] pupil,
                int device )

#################################################
# P-Class Target
#################################################
cdef class Target:
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


cdef class Sensors:
    cdef sutra_sensors *sensors

    cdef sensors_initgs(self,np.ndarray[dtype=np.float32_t] xpos,
                             np.ndarray[dtype=np.float32_t] ypos,
                             np.ndarray[dtype=np.float32_t] Lambda,
                             np.ndarray[dtype=np.float32_t] mag,
                             np.ndarray[dtype=np.int64_t  ] size,
                             np.ndarray[dtype=np.float32_t] noise,
                             np.ndarray[dtype=np.int64_t  ] seed)

    #cdef sensors_initarr(self, int n, int *phasemap, float *offset, float *pupil,
     #       int *isvalid, int *validsubx, int *validsuby, float* flux_foc, 
     #       int *istart, int *jstart, float *ftkern_sinc, int *hrmap, 
     #       int *binmap, float *halfxy)

    cdef sensors_initarr(self,int n, Param_wfs wfs, Param_geom geom)
    cdef sensors_addlayer(self,int i, bytes type, float alt, float xoff, float yoff)
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

#################################################
# P-Class (parametres) Param_wfs
#################################################
cdef class Param_wfs:
    cdef readonly str type
    """type of wfs : "sh" or "pyr"."""
    cdef readonly long nxsub
    """linear number of subaps."""
    cdef readonly long  npix
    """number of pixels per subap."""
    cdef readonly float pixsize
    """pixel size (in arcsec) for a subap."""
    cdef readonly float Lambda
    """observation wavelength (in µm) for a subap."""
    cdef readonly float optthroughput
    """wfs global throughput."""
    cdef readonly float fracsub
    """minimal illumination fraction for valid subaps."""
    cdef readonly long openloop
    """1 if in "open-loop" mode (i.e. does not see dm)."""
    cdef readonly float fssize
    """size of field stop in arcsec."""
    cdef readonly str fstop
    """size of field stop in arcsec."""

    cdef readonly int atmos_seen
    """1 if the WFS sees the atmosphere layers"""
    #TODO pointer dms_seen;      // index of dms seen by the WFS


#target kwrd
    cdef readonly float xpos
    """guide star x position on sky (in arcsec)."""
    cdef readonly float ypos
    """guide star x position on sky (in arcsec)."""
    cdef readonly float gsalt
    """altitude of guide star (in m) 0 if ngs."""
    cdef readonly float gsmag
    """magnitude of guide star."""
    cdef readonly float zerop
    """detector zero point."""
    cdef readonly float noise
    """desired noise : < 0 = no noise / 0 = photon only / > 0 photon + ron."""

    cdef readonly float kernel      # 
    cdef np.ndarray _ftkernel #(float*)

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

#internal kwrd
    cdef long  _pdiam          # pupil diam for a subap (in pixel)
    cdef long  _Nfft           # array size for fft for a subap (in pixel)
    cdef long  _Ntot           # total size of hr image for a subap (in pixel)
    cdef long  _nrebin         # rebin factor from hr to binned image for a subap 
    cdef long  _nvalid         # number of valid subaps

    cdef float   _nphotons     # number of photons per subap
    cdef float   _subapd       # subap diameter (m)
    cdef readonly np.ndarray _fluxPerSub   # fraction of nphotons per subap
    cdef float   _qpixsize     # quantum pixel size for the simulation

    cdef readonly np.ndarray _istart    # (int*) x start indexes for cutting phase screens 
    cdef readonly np.ndarray _jstart    # (int*) y start indexes for cutting phase screens 
    #cdef np.ndarray _validsubs    # (i,j) indices of valid subaps
    cdef readonly np.ndarray _validsubsx    #(int*) indices of valid subaps along axis x
    cdef readonly np.ndarray _validsubsy    #(int*) indices of valid subaps along axis y
    cdef readonly np.ndarray _isvalid      #(int*) array of 0/1 for valid subaps
    cdef readonly np.ndarray _phasemap     #(int*) array of pixels transform from phase screen into
                         # subaps phase screens
    cdef readonly np.ndarray _hrmap        #(int*) array of pixels transform from minimal FoV image to (in case type is "sh" or "geo"
    cdef np.ndarray _sincar        #(float*) array of pixels transform from minimal FoV image to (in case type is "pyr" or "roof")
                         # full FoV image (array of 0 if the same)
    cdef readonly np.ndarray _binmap       #(int*) array of pixels transform from full FoV hr images to
                         # binned images
    cdef readonly np.ndarray _halfxy       #(float*) phase offset for 1/2 pixel shift in (x,y)

    cdef readonly np.ndarray _submask       #(float*) fieldstop for each subap

    cdef float* _lgskern      # lgs kernels for each subap
    cdef float* _profna       # sodium profile
    cdef float* _altna        # corresponding altitude
    cdef float* _prof1d       # hr profile
    cdef float* _profcum      # hr profile cumulated
    cdef float* _beam         # 1d beam function
    cdef float* _ftbeam       # 1d beam function fft
    cdef float* _azimuth      # angles of rotation for each spot

# pyramid-nly kwrds
    cdef readonly float pyr_ampl
    """pyramid wfs modulation amplitude radius [arcsec]."""
    cdef readonly long pyr_npts
    """total number of point along modulation circle [unitless]."""
    cdef np.ndarray pyr_pos    # positions for modulation, overwrites ampl and npts [arcsec]
    cdef readonly str pyr_loc
    """Location of modulation, before/after the field stop.
    valid value are "before" or "after" (default "after")."""
    cdef readonly str pyrtype
    """Type of pyramid, either 0 for "Pyramid" or 1 for "RoofPrism"."""

# pyramid internal kwrds
    cdef readonly np.ndarray _pyr_offsets   #(float*)
    cdef readonly np.ndarray _pyr_cx   #(int*)
    cdef readonly np.ndarray _pyr_cy   #(int*)


#################################################
# P-Class (parametres) Param_dm
#################################################
cdef class Param_dm:
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


#################################################
# P-Class (parametres) Param_loop
#################################################
cdef class Param_loop:
    cdef long niter     # number of iterations
    cdef float ittime   # iteration time (in sec)


