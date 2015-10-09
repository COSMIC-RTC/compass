include "../par.pxi"

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

from cpython.string cimport PyString_AsString

IF USE_MPI == 1:
    print "using mpi"
    #from mpi4py import MPI
    from mpi4py cimport MPI
    #cimport mpi4py.MPI as MPI
    # C-level cdef, typed, Python objects
    from mpi4py cimport mpi_c as mpi
    #from mpi4py cimport libmpi as mpi
    #cimport mpi4py.libmpi as mpi


from libc.math cimport sin

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
        #carma_obj[float] *d_pupil
        carma_context *current_context

        sutra_atmos(carma_context *context, int nscreens, float *r0, long *size,
                    long *size2, float *altitude, float *windspeed, float *winddir,
                    float *deltax, float *deltay, int device)
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
        sutra_lgs *d_lgs #   the lgs object
        carma_obj[float] *object #   the object intensity map
        carma_obj[float] *d_pupil #   the pupil mask
        carma_obj[float] *d_image #   the resulting image for this target
        carma_obj[float] *d_leimage #   the long exposure image for this target
        carma_obj[cuFloatComplex] *d_amplipup #   the complex amplitude in the pupil plane
        carma_obj[float] *d_phasepts #   the valid phase points in the pupil (target only)
        carma_obj[int] *d_wherephase #   positions of valid phase points in the pupil (target only)
        #type_screen unknown cf sutra_DM
        map[pair[string,int], float] xoff #   x reference for raytracing
        map[pair[string,int], float] yoff #   y reference for raytracing
        carma_context *current_context
        
        int add_layer(char* l_type, float alt, float xoff, float yoff)
        int init_strehlmeter()
        int raytrace(sutra_atmos *atmos)
        int raytrace(sutra_dms *ydms, int rst)
        int comp_image(int puponly)
        int comp_strehl()
        
    cdef cppclass sutra_target:
        int ntargets
        vector.vector[sutra_source *] d_targets

        sutra_target(carma_context *context, int ntargets, float *xpos, float *ypos,
      float *Lambda, float *mag, float zerop, long *sizes, float *pupil, int Npts, int device)

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
        int device
        carma_context *current_context
        size_t nsensors() 
        
        vector[sutra_wfs *] d_wfs
        map[vector[int],cufftHandle*] campli_plans
        map[vector[int],cufftHandle*] fttotim_plans
        map[vector[int],cufftHandle*] ftlgskern_plans

        carma_obj[cuFloatComplex] *d_camplipup
        carma_obj[cuFloatComplex] *d_camplifoc
        carma_obj[cuFloatComplex] *d_fttotim
        carma_obj[cuFloatComplex] *d_ftlgskern
        carma_obj[float] *d_lgskern

        sutra_sensors(carma_context *context, const char** type, int nwfs, long *nxsub,
          long *nvalid, long *npix, long *nphase, long *nrebin, long *nfft,
          long *ntot, long *npup, float *pdiam, float *nphot, int *lgs, int device)
        sutra_sensors(carma_context *context, int nwfs, long *nxsub, long *nvalid,
          long *nphase, long npup, float *pdiam, int device)

        int sensors_initgs(float *xpos, float *ypos, float *Lambda, float *mag, float zerop,
                   long *size, float *noise, long *seed)
        int sensors_initgs(float *xpos, float *ypos, float *Lambda, float *mag,float zerop,
                   long *size, float *noise)
        int sensors_initgs(float *xpos, float *ypos, float *Lambda, float *mag,float zerop,
                   long *size)
        int allocate_buffers()
        int define_mpi_rank(int rank, int size)


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

        int rank
        int offset
        int nvalid_tot
        int* displ_bincube
        int* count_bincube

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
          float *pupil, float *fluxPerSub, int *validsubsx,
          int *validsubsy, int *istart, int *jstart, cuFloatComplex *kernel)
        int wfs_initarrays(int *phasemap, float *offsets, float *pupil, float *fluxPerSub,
          int *validsubsx, int *validsubsy, int *istart, int *jstart)
        int wfs_initarrays(cuFloatComplex *halfxy, cuFloatComplex *offsets,
          float *focmask, float *pupil, int *cx, int *cy,
          float *sincar, int *phasemap, int *validsubsx, int *validsubsy)
        int wfs_initgs(sutra_sensors *sensors, float xpos, float ypos, float Lambda, float mag, long size,
          float noise, long seed)
        int load_kernels(float *lgskern)
        int sensor_trace(sutra_atmos *yatmos)
        int sensor_trace(sutra_atmos *atmos, sutra_dms *ydms)
        int sensor_trace(sutra_dms *ydm, int rst)
        int comp_image()
        int comp_generic()
        int define_mpi_rank(int rank, int size)
        int allocate_buffers(sutra_sensors *sensors)


#################################################
# C-Class sutra_wfs_sh
#################################################
cdef extern from "sutra_wfs_sh.h":
    cdef cppclass sutra_wfs_sh(sutra_wfs):
        int  wfs_initarrays(int *phasemap, int *hrmap, int *binmap, float *offsets,
            float *pupil, float *fluxPerSub, int *validsubsx,
            int *validsubsy, int *istart, int *jstart, cuFloatComplex *kernel)
        int fill_binimage(int async)
        int slopes_geom(int type, float *slopes)
        int slopes_geom(int type)


#################################################
# C-Class sutra_wfs_pyr_roof
#################################################
cdef extern from "sutra_wfs_pyr.h":
    cdef cppclass sutra_wfs_pyr(sutra_wfs):
        int wfs_initarrays(cuFloatComplex *halfxy, cuFloatComplex *offsets,
      float *focmask, float *pupil, int *cx, int *cy,
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
            float *fluxPerSub, int *validsubsx, int *validsubsy)


#################################################
# C-Class sutra_centroider
#################################################
cdef extern from "sutra_centroider.h":
    cdef cppclass sutra_centroider:
        int device
        sutra_wfs *wfs
        int nwfs
        int nvalid

        float offset
        float scale

        carma_context *current_context

        bool is_type(string typec)

        string get_type()

        int get_cog(carma_streams *streams, float *cube, float *subsum, float *centroids,
                int nvalid, int npix, int ntot)
        int get_cog(float *subsum, float *slopes)
        int get_cog()



    int convert_centro(float *d_odata, float *d_idata, float offset, float scale,
    int N, carma_device *device)


#################################################
# C-Class sutra_centroider_tcog
#################################################
cdef extern from "sutra_centroider_tcog.h":
    cdef cppclass sutra_centroider_tcog(sutra_centroider):
        int set_threshold(float threshold)
        bool is_type(string typec)
        #string get_type() 

#################################################
# C-Class sutra_centroider_corr
#################################################
cdef extern from "sutra_centroider_corr.h":
    cdef cppclass sutra_centroider_corr(sutra_centroider):
        int init_bincube()
        bool is_type(string typec)
        string get_type() 

        int init_corr(int isizex, int isizey, float *interpmat)
        int load_corr(float *corr, float *corr_norm, int ndim)


#################################################
# C-Class sutra_centroider_wcog
#################################################
cdef extern from "sutra_centroider_wcog.h":
    cdef cppclass sutra_centroider_wcog(sutra_centroider):
        int npix
        carma_obj[float] *d_weights

        sutra_centroider_wcog(carma_context *context, sutra_sensors *sensors,
            int nwfs, long nvalid, float offset, float scale, int device)
        sutra_centroider_wcog(const sutra_centroider& centroider)
        string get_type()
        int init_weights()
        int load_weights(float *weights, int ndim)
        int get_cog(carma_streams *streams, float *cube, float *subsum, 
            float *centroids, int nvalid, int npix, int ntot)
        int get_cog(float *subsum, float *slopes)
        int get_cog()



#################################################
# Dynamic casts
#################################################
cdef extern from *:
    sutra_wfs_geom* dynamic_cast_wfs_geom_ptr "dynamic_cast<sutra_wfs_geom*>" (sutra_wfs*) except NULL
    sutra_wfs_sh* dynamic_cast_wfs_sh_ptr "dynamic_cast<sutra_wfs_sh*>" (sutra_wfs*) except NULL
    sutra_wfs_pyr_pyr4* dynamic_cast_wfs_pyr_pyr4_ptr "dynamic_cast<sutra_wfs_pyr_pyr4*>" (sutra_wfs*) except NULL
    sutra_wfs_pyr_roof* dynamic_cast_wfs_pyr_roof_ptr "dynamic_cast<sutra_wfs_pyr_roof*>" (sutra_wfs*) except NULL
    sutra_centroider_tcog* dynamic_cast_centroider_tcog_ptr "dynamic_cast<sutra_centroider_tcog*>" (sutra_centroider*) except NULL
    sutra_centroider_corr* dynamic_cast_centroider_corr_ptr "dynamic_cast<sutra_centroider_corr*>" (sutra_centroider*) except NULL
    sutra_centroider_wcog* dynamic_cast_centroider_wcog_ptr "dynamic_cast<sutra_centroider_wcog*>" (sutra_centroider*) except NULL
    sutra_controller_generic* dynamic_cast_controller_generic_ptr "dynamic_cast<sutra_controller_generic*>" (sutra_controller*) except NULL
    sutra_controller_geo* dynamic_cast_controller_geo_ptr "dynamic_cast<sutra_controller_geo*>" (sutra_controller*) except NULL
    sutra_controller_ls* dynamic_cast_controller_ls_ptr "dynamic_cast<sutra_controller_ls*>" (sutra_controller*) except NULL
    sutra_controller_mv* dynamic_cast_controller_mv_ptr "dynamic_cast<sutra_controller_mv*>" (sutra_controller*) except NULL
    sutra_controller_cured* dynamic_cast_controller_cured_ptr "dynamic_cast<sutra_controller_cured*>" (sutra_controller*) except NULL
    sutra_controller_kalman* dynamic_cast_controller_kl_ptr "dynamic_cast<sutra_controller_kalman*>" (sutra_controller*) except NULL


cdef extern from "sutra_dm.h":
#################################################
# C-Class sutra_dms
#################################################
    ctypedef pair[string, float] type_screen

    cdef cppclass sutra_dms:
        vector[sutra_dm *] d_dms
        vector[type_screen] d_type

        sutra_dms(int ndm)

        int add_dm(carma_context *context, const char* type, float alt, long dim,
          long ninflu, long influsize, long ninflupos, long n_npoints,
          float push4imat, int device)
        int remove_dm(string type, float alt)

        int get_inddm(string type, float alt)

        int ndm()
        int nact_total()

#################################################
# C-Class sutra_dm
#################################################
    cdef cppclass sutra_dm:
        int device
        string type
        long ninflu
        long influsize
        long dim
        float push4imat

        sutra_phase *d_shape

        carma_obj[float] *d_comm

        carma_obj[float] *d_influ # if relevant

        # pzt 
        carma_obj[int] *d_influpos
        carma_obj[int] *d_npoints
        carma_obj[int] *d_istart
        carma_obj[int] *d_xoff
        carma_obj[int] *d_yoff
        carma_obj[int] *d_pos # Convolution preprocess
        carma_obj[float] *d_KLbasis
        carma_obj[float] *d_kernconv # Convolution preprocess
        carma_obj[float] *d_mapactu # Convolution process
        carma_obj[cuFloatComplex] *d_ftkernconv # Convolution process
        carma_obj[cuFloatComplex] *d_ftmapactu # Convolution process

        #zernike
        carma_obj[float] *d_coeffs
        carma_obj[float] *d_mask
        carma_obj[float] *d_zr
        carma_obj[float] *d_ztheta

        #sutra_kl *d_kl

        carma_context *current_context




        int pzt_loadarrays(float *influ, float *influ2, tuple_t[float] *influ3,
                            int *influpos, int *influpos2, int *npoints, int *istart,
                            int *xoff, int *yoff, float *kernconv)
        int kl_loadarrays(float *rabas, float *azbas, int *ord, float *cr, float *cp)

        int reset_shape()

        int comp_oneactu(int nactu, float ampli)

        int comp_shape()
        int comp_shape(float *comm)
        int compute_KLbasis(float *xpos, float *ypos, int *indx, long dim,
            float norm, float ampli)

#################################################
# C-Class sutra_controller
#################################################
cdef extern from "sutra_controller.h":

    cdef cppclass sutra_controller:
#allocation of d_centroids and d_com
        sutra_controller(carma_context* context, int nslope, int nactu, float delay,
            sutra_dms *dms, char **type, float *alt, int ndm)

        string get_type()

#!!!! YOU MUST set d_centroids before calling it!!!!
        int comp_com()

#It is better to have something like this (+protected d_centroids):
#virtual int comp_com (carma_obj<float> *new_centroids)=0
#it would imply copy, but would be much safer

        inline int nactu()
        inline int nslope()


        #TODO cublasHandle_t cublas_handle()

        int set_perturbcom(float *perturb, int N)
        int set_openloop(int open_loop_status)
        int comp_voltage()
        int syevd_f(char meth, carma_obj[float] *d_U, carma_host_obj[float] *h_eingenvals)
        int invgen(carma_obj[float] *d_mat, float cond, int job)
        int command_delay()

        #I would propose to make them protected (+ proper
        #set of fuctions). It could make life easier!
        #But we should discuss it
        int cpt_pertu
        int open_loop
        float delay
        float a # Coefficient for linear interpolation on command buffer to allow non-integer delay
        float b # Coefficient for linear interpolation on command buffer to allow non-integer delay
        float c # Coefficient for linear interpolation on command buffer to allow non-integer delay
        vector[sutra_dm *] d_dmseen
        carma_obj[float] *d_centroids # current centroids
        carma_obj[float] *d_com # current command
        carma_obj[float] *d_perturb # perturbation command buffer
        carma_obj[float] *d_voltage # commands sent to mirror
        carma_obj[float] *d_com1 # commands k-1
        carma_obj[float] *d_com2 # commands k-2

        carma_streams *streams

        int device
        carma_context *current_context



#################################################
# C-Class sutra_controller_generic
#################################################
cdef extern from "sutra_controller_generic.h":
    cdef cppclass sutra_controller_generic:
        carma_obj[float] *d_matE
        carma_obj[float] *d_cmat
        carma_obj[float] *d_gain
        carma_obj[float] *d_decayFactor
        carma_obj[float] *d_compbuff
        string command_law

        sutra_controller_generic(carma_context *context, long nvalid, long nactu,
          float delay, sutra_dms *dms, char **type, float *alt, int ndm)
        sutra_controller_generic(const sutra_controller_generic& controller)

        string get_type()
        string get_commandlaw()
        int set_decayFactor(float *decayFactor)
        int set_mgain(float *gain)
        int set_cmat(float *cmat)
        int set_matE(float *matE)
        int set_commandlaw(string law)
        int comp_com()


#################################################
# C-Class sutra_controller_geo
#################################################
cdef extern from "sutra_controller_geo.h":
    cdef cppclass sutra_controller_geo:
        float gain
        long Nphi

        carma_obj[float] *d_gain
        carma_obj[float] *d_proj
        carma_obj[double] *d_phi
        carma_obj[int] *d_indx_pup
        # TODO carma_sparse_obj[double] *d_IFsparse
        carma_obj[float] *d_geocov
        carma_obj[double] *d_compdouble
        carma_obj[float] *d_compfloat

        sutra_controller_geo(carma_context *context, long nactu, long Nphi,
            float delay, sutra_dms *dms, char **type, float *alt, int ndm)
        sutra_controller_geo(const sutra_controller_geo& controller)

        string get_type()

        int set_gain(float gain)
        int load_mgain(float *mgain)
        int comp_dphi(sutra_source *target)
        int comp_com()
        int init_proj(sutra_dms *dms, int *indx_dm, float *unitpervolt, int *indx_pup)
        int init_proj_sparse(sutra_dms *dms, int *indx_dm, float *unitpervolt, int *indx_pup)



#################################################
# C-Class sutra_controller_ls
#################################################
cdef extern from "sutra_controller_ls.h":
    cdef cppclass sutra_controller_ls:
        float gain

        carma_obj[float] *d_imat
        carma_obj[float] *d_cmat
        carma_obj[float] *d_gain

# svd computations
        carma_obj[float] *d_eigenvals
        carma_host_obj[float] *h_eigenvals
        carma_obj[float] *d_U

# loop components
        carma_obj[float] *d_cenbuff # centroids circular buffer
        carma_obj[float] *d_err # current error

# Modal optimization components
        int is_modopti # Flag for using modal optimization
        int nrec # Number of recorded open slopes measurements
        int nmodes #Number of modes
        float gmin # Gain min
        float gmax # Gain max
        int ngain # Number of tested gains between gmin and gmax
        float Fs # Sampling frequency
        int cpt_rec # Counter for modal gains refresh
        carma_obj[float] *d_M2V # Modes to Volts matrix
        carma_obj[float] *d_S2M # Slopes to Modes matrix
        carma_obj[float] *d_slpol # Open-loop measurements buffer, recorded and loaded from Yorick
        carma_obj[float] *d_Hcor # Transfer function
        carma_obj[float] *d_com1 # Command k-1 for POLC
        carma_obj[float] *d_com2 # Command k-2 for POLC
        carma_obj[float] *d_compbuff # Buffer for POLC computation
        carma_obj[float] *d_compbuff2 # Buffer for POLC computation

        sutra_controller_ls(carma_context *context, long nvalid, long nactu,
          float delay, sutra_dms *dms, char **type, float *alt, int ndm)
        sutra_controller_ls(const sutra_controller_ls& controller)

        string get_type()

        int svdec_imat()
        int build_cmat(int nfilt, bool filt_tt)
        int build_cmat(int nfilt)
        int build_cmat_modopti()
        int frame_delay()
        int comp_com()
        int set_gain(float gain)
        int set_mgain(float *mgain)
        int set_cmat(float *cmat)
        int set_delay(float delay)
        int init_modalOpti(int nmodes, int nrec, float *M2V, float gmin, float gmax, int ngain, float Fs)
        int loadOpenLoopSlp(float *ol_slopes)
        int modalControlOptimization()
        int compute_Hcor() 


#################################################
# C-Class sutra_controller_mv
#################################################
cdef extern from "sutra_controller_mv.h":
    cdef cppclass sutra_controller_mv:

        float gain

        carma_obj[float] *d_imat
        carma_obj[float] *d_cmat
        carma_obj[float] *d_gain

# Cphim & Cmm features
        carma_obj[float] *d_covmat
        carma_obj[float] *d_KLbasis
        carma_obj[float] *d_noisemat
        carma_obj[float] *d_Cmm
        carma_obj[float] *d_Cphim
# svd computations
        carma_host_obj[float] *h_Cmmeigenvals
        carma_host_obj[float] *h_eigenvals
#carma_obj[float] *d_U

# loop components
        carma_obj[float] *d_cenbuff # centroids circular buffer
        carma_obj[float] *d_com1 # commands k-1 (for POLC)
        carma_obj[float] *d_com2 # commands k-2 (for POLC)
        carma_obj[float] *d_compbuff # Buffer for computations
        carma_obj[float] *d_compbuff2
        carma_obj[float] *d_olmeas # Open-loop measurements for POLC
        carma_obj[float] *d_err # current error

        #TODO cublasHandle_t cublas_handle

        sutra_controller_mv(carma_context *context, long nvalid, long nactu,
          float delay, sutra_dms *dms, char **type, float *alt, int ndm)
        sutra_controller_mv(const sutra_controller_mv& controller)

        string get_type()

        int svdec_imat()
        int build_cmat(const char *dmtype, char *method)
        int build_cmat(float cond)
        int frame_delay()
        int comp_com()
        int set_gain(float gain)
        int set_mgain(float *mgain)
        int set_delay(float delay)
        int set_cmat(float *cmat)
        # Florian features
        int load_noisemat(float *noise)
        int do_covmat(sutra_dm *ydm,char *method, int *indx_pup, long dim, float *xpos, float *ypos, long Nkl, float norm, float ampli)
        int do_geomat(carma_obj[float] *d_geocov, carma_obj[float] *d_IF, long n_pts, float ampli)
        int piston_filt(carma_obj[float] *d_statcov)
        int piston_filt_cphim(carma_obj[float] *d_cphim, float *F)
        int filter_cphim(float *F, float *Nact)
        int filter_cmat(float cond)
        int invgen(carma_obj[float] *d_mat, float cond, int job)
        int invgen(carma_obj[float] *d_mat, carma_host_obj[float] *h_eigen, float cond)
        int invgen_cpu(carma_obj[float] *d_mat, carma_host_obj[float] *h_eigen, float cond)
        # int do_statmat(float *statcov,long dim, float *xpos, float *ypos, float norm, carma_device *device)
        int DDiago(carma_obj[float] *d_statcov, carma_obj[float] *d_geocov)
        int load_covmat(float *covmat)
        int load_klbasis(float *klbasis)
        int compute_Cmm(sutra_atmos *atmos, sutra_sensors *sensors, double *L0, double *cn2, double *alphaX, double *alphaY, double diamTel, double cobs)
        int compute_Cphim(sutra_atmos *atmos,
                sutra_sensors *sensors, sutra_dms *dms, double *L0, double *cn2,
                double *alphaX, double *alphaY, double *X, double *Y, double *xactu,
                double *yactu, double diamTel, double *k2, long *NlayerDm,
                long *indLayerDm, double FoV, double *pitch, double *alt_dm)


#################################################
# C-Class sutra_controller_cured
#################################################
cdef extern from "sutra_controller_cured.h":
    cdef cppclass sutra_controller_cured:
        float gain
        int   ndivs #number of subdivision levels for cured
        bool  tt_flag #flag for separate tt

        # data for CuReD */
        carma_host_obj[float] *h_centroids
        carma_host_obj[float] *h_err
        carma_obj[float] *d_err # current error
        carma_obj[float] *d_cenbuff # centroids circular buffer

        # data for CuReD */
        carma_obj[float] *d_imat

        # structures needed to run CuReD */
        #sysCure* h_syscure
        void* h_syscure
        #parCure* h_parcure
        void* h_parcure

        sutra_controller_cured(carma_context *context, long nvalid, long nactu,
                 float delay, sutra_dms *dms, char **type, float *alt, int ndm)
        sutra_controller_cured(const sutra_controller_cured& controller)

        string get_type()
        int set_gain(float gain)

        int comp_com()

        int init_cured(int nxsubs, int *isvalid, int ndivs, int tt)
        int frame_delay()
        int set_delay(float delay)


#################################################
# C-Class sutra_controller_kalman
#################################################
cdef extern from "sutra_controller_kalman.h":
    cdef cppclass sutra_controller_kalman:
        sutra_controller_kalman(carma_context* context, int nvalid_, int nactu_, sutra_dms *dms, char **type, float *alt, int ndm)

        void init_kalman(carma_host_obj[float]& D_Mo, carma_host_obj[float]& N_Act, carma_host_obj[float]& PROJ, bool is_zonal, bool is_sparse, bool is_GPU)

        double gettime()
        double gettime_op1()
        double gettime_op2()
        double gettime_op3()

        void calculate_gain(float bruit, carma_host_obj[float]& SigmaV,
                 carma_host_obj[float]& atur, carma_host_obj[float]& btur)

        int set_gain(float k_W)

        #TODO cusparseHandle_t cusparseHandle
        #kp_kalman_core_sparse* core_sparse
        #kp_kalman_core_full* core_full
        bool isGPU
        bool isSparse
        bool isZonal
        bool isInit
        bool isGainSet
        float gain



#################################################
# C-Class sutra_rtc
#################################################
cdef extern from "sutra_rtc.h":
    cdef cppclass sutra_rtc:
        int device

        vector[sutra_centroider*] d_centro
        vector[sutra_controller*] d_control

        carma_context *current_context

        sutra_rtc(carma_context *context)
        #sutra_rtc(const sutra_rtc& yrtc)

        int add_centroider(sutra_sensors *sensors, int nwfs, long nvalid, float offset, float scale, long device,
             char *typec)
        int rm_centroider()
        int add_controller_geo(int nactu, int Nphi, float delay, long device,
                sutra_dms *dms, char **type_dmseen, float *alt, int ndm)
        int add_controller(int nactu, float delay, long device, const char *typec,
                sutra_dms *dms, char **type_dmseen, float *alt, int ndm)
        int rm_controller()

        int do_imat(int ncntrl, sutra_dms *ydms)
        #int do_imatkl(int ncntrl, sutra_dms *ydms)
        #int do_imatkl4pzt(int ncntrl, sutra_dms *ydms)
        int do_imat_geom(int ncntrl, sutra_dms *ydm, int type)

        int do_centroids()
        int do_centroids(int ncntrl)
        int do_centroids(int ncntrl, bool imat)
        int do_centroids_geom(int ncntrl)
        int do_control(int ncntrl)
        int apply_control(int ncntrl, sutra_dms *ydm)

#################################################
# C-Class sutra_lgs
#################################################
cdef extern from "sutra_lgs.h":
    cdef cppclass sutra_lgs:
        int device
        long nvalid
        long npix
        long nmaxhr
        float hg
        float h0
        float deltah
        float pixsize
        long nprof

        cufftHandle *ftlgskern_plan
        carma_obj[float] *d_doffaxis
        carma_obj[float] *d_azimuth
        carma_obj[float] *d_prof1d
        carma_obj[float] *d_profcum
        carma_obj[cuFloatComplex] *d_prof2d
        carma_obj[float] *d_beam
        carma_obj[cuFloatComplex] *d_ftbeam
        carma_obj[float] *d_lgskern
        carma_obj[cuFloatComplex] *d_ftlgskern

        carma_context *current_context

        sutra_lgs(carma_context *context, sutra_sensors *sensors, long nvalid, long npix, long nmaxhr)
        sutra_lgs(const sutra_lgs& lgs)

        int lgs_init(int nprof, float hg, float h0, float deltah, float pixsie,
            float *doffaxis, float *prof1d, float *profcum, float *beam,
            cuFloatComplex *ftbeam, float* azimuth)
        int load_prof(float *prof1d, float *profcum, float hg, float h0, float deltah)
        int lgs_update(carma_device *device)
        int lgs_makespot(carma_device *device, int nin)
        int load_kernels(float *lgskern, carma_device *device)




#################################################
# P-Class Atmos
#################################################
cdef class Atmos:
    cdef sutra_atmos *s_a
    cdef chakra_context context
    cdef realinit(self,chakra_context ctxt,int nscreens,
                np.ndarray[ndim=1,dtype=np.float32_t] r0,
                np.ndarray[ndim=1,dtype=np.int64_t] size,
                np.ndarray[ndim=1,dtype=np.float32_t] altitude,
                np.ndarray[ndim=1,dtype=np.float32_t] windspeed,
                np.ndarray[ndim=1,dtype=np.float32_t] winddir,
                np.ndarray[ndim=1,dtype=np.float32_t] deltax,
                np.ndarray[ndim=1,dtype=np.float32_t] deltay,
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

    cpdef dmtrace(self,int ntar, Dms dms, int reset=?)

#################################################
# P-Class Sensors
#################################################
cdef class Sensors:
    
    cdef sutra_sensors *sensors
    cdef sensors_initgs(self,np.ndarray[ndim=1,dtype=np.float32_t] xpos,
                             np.ndarray[ndim=1,dtype=np.float32_t] ypos,
                             np.ndarray[ndim=1,dtype=np.float32_t] Lambda,
                             np.ndarray[ndim=1,dtype=np.float32_t] mag,
                             float zerop,
                             np.ndarray[ndim=1,dtype=np.int64_t  ] size,
                             np.ndarray[ndim=1,dtype=np.float32_t] noise,
                             np.ndarray[ndim=1,dtype=np.int64_t  ] seed)


    cdef sensors_initarr(self,int n, Param_wfs wfs, Param_geom geom)
    cdef sensors_addlayer(self,int i, bytes type, float alt, float xoff, float yoff)
    cdef _get_bincube(self, int n)
    cdef _get_binimg(self, int n)
    cdef _get_slopesDims(self,int n)
    cdef _get_slopes(self, int n)
    cdef slopes_geom(self,int nsensors, int t)
    cpdef sensors_trace(self,int n, str type_trace, Atmos atmos=?, Dms dms=?, int rst=?)
    IF USE_MPI==1:
        cpdef gather_bincube(self,int n)
        cpdef gather_bincube_cuda_aware(self,int n)
        cpdef Bcast_dscreen(self)
        cpdef Bcast_dscreen_cuda_aware(self)
    cdef _get_rank(self,int n)

    #for profiling purpose
    '''
    cdef gather_bincube_prof(self,int n)
    cdef wait1_prof(self)
    cdef wait2_prof(self)
    cdef d2h_prof(self,float *ptr,n)
    cdef h2d_prof(self,float *ptr,n)
    cdef gather_prof(self,float *ptr, int size, int *count, int *disp)
    '''

    cdef  _get_hrmap(self, int n)
    #cdef getDims(self)

#################################################
# P-Class Rtc
#################################################
cdef class Rtc:
    cdef sutra_rtc *rtc
    cdef int device
    cdef int use_brama
    #cdef sensors_initbcube(self,int ncentro)
    cpdef getcentroids(self,int ncontrol, Sensors g_wfs=?, int nwfs=?)
    cpdef docentroids(self,int ncontrol=?)
    cdef init_proj(self,int i,Dms dms,np.ndarray[ndim=1,dtype=np.int32_t] indx_dm,
            np.ndarray[ndim=1,dtype=np.float32_t] unitpervolt, 
            np.ndarray[ndim=1,dtype=np.int32_t] indx_pup)
    cdef init_modalOpti(self,int ncontro,int nmodes,int nrec, np.ndarray[ndim=2,dtype=np.float32_t] M2V,
            float gmin, float gmax, int ngain, float Fs)
    cdef loadOpenLoop(self,int ncontro, np.ndarray[ndim=2, dtype=np.float32_t] ol_slopes)
    cdef modalControlOptimization(self,int ncontro)
    cdef set_gain(self,int ncontro, float gain)
    cdef set_mgain(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] mgain)
    cpdef set_imat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data)
    cpdef get_imat(self, int ncontro)
    cpdef set_cmat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data)
    cpdef get_cmat(self,int ncontro)
    cpdef get_cphim(self, int ncontro)
    cpdef get_cmm(self,int ncontro)
    cdef set_decayFactor(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] decay)
    cdef set_matE(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] matE)


    cdef doimat_geom(self, int ncontro, Dms g_dms,int geom)
    cdef doimat(self, int ncontro, Dms g_dms)
    cdef sensors_compslopes(self, int ncentro, int nmax=?, float thresh=?)
    cdef add_controller(self, int nactu, float delay, bytes type_control, Dms dms,
                 char **type_dmseen, np.ndarray[ndim=1,dtype=np.float32_t] alt,
                 int ndm, long Nphi=?)


    cdef imat_svd(self,int ncontro)
    cdef setU(self,int ncontro,np.ndarray[ndim=2,dtype=np.float32_t] U)
    cdef setEigenvals(self, int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] eigenvals)
    cdef getU(self, int ncontro)
    cdef getEigenvals(self,int ncontro)
    cpdef getCenbuff(self, int ncontro)
    cdef getErr(self,int ncontro)
    cpdef getCom(self,int ncontro)
    cpdef getolmeas(self,int ncontro)
    cpdef getVoltage(self,int ncontro)
    cpdef getCentroids(self,int ncontro)
    cdef buildcmat(self,int ncontro,int nfilt, int filt_tt=?)
    cdef buildcmatmv(self,int ncontro,float cond)
    cdef loadnoisemat(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] N)
    cpdef docontrol(self,int ncontro)
    cpdef applycontrol(self,int ncontro,Dms dms)


#################################################
# P-Class Dms
#################################################
cdef class Dms:
    cdef sutra_dms *dms
    cdef int device
    cdef add_dm(self, bytes type_dm, float alt, long dim, long ninflu, 
                long influsize, long ninflupos, long npts, float puhs4imat,
                int device=?)
    cdef remove_dm(self,bytes type_dm, float alt)

    cdef load_pzt(self, float alt,
                    np.ndarray[ndim=3,dtype=np.float32_t] influ,
                    np.ndarray[ndim=1,dtype=np.int32_t] influpos,
                    np.ndarray[ndim=1,dtype=np.int32_t] npoints,
                    np.ndarray[ndim=1,dtype=np.int32_t] istart,
                    np.ndarray[ndim=1,dtype=np.int32_t] xoff,
                    np.ndarray[ndim=1,dtype=np.int32_t] yoff,
                    np.ndarray[ndim=2,dtype=np.float32_t] kern)

    #TODO dims of arrays
    cdef load_kl(self,float alt, np.ndarray[ndim=1,dtype=np.float32_t] rabas,
                    np.ndarray[ndim=1,dtype=np.float32_t] azbas,
                    np.ndarray[ndim=1,dtype=np.int32_t] ord,
                    np.ndarray[ndim=1,dtype=np.float32_t] cr,
                    np.ndarray[ndim=1,dtype=np.float32_t] cp)


    cdef load_tt(self,float alt, np.ndarray[ndim=3,dtype=np.float32_t] influ)

    cdef set_comm(self,bytes type_dm,float alt,
                    np.ndarray[ndim=1,dtype=np.float32_t] comm)
    cdef shape_dm(self,bytes type_dm,float alt)

    cdef computeKLbasis(self, bytes type_dm, float alt, 
        np.ndarray[ndim=1,dtype=np.float32_t] xpos, np.ndarray[ndim=1,dtype=np.float32_t] ypos,
        np.ndarray[ndim=1,dtype=np.int32_t] indx_pup, long dim, float norm, float ampli)
    cdef get_KLbasis(self,bytes type_dm, float alt)

    cpdef getComm(self,bytes type_dm,float alt)
    cpdef getInflu(self,bytes type_dm,float alt)
    cpdef comp_oneactu(self,bytes type_dm, float alt, int nactu, float ampli)

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
    cdef readonly np.ndarray _mpupil   # medium pupil (part of the guard band)
    cdef np.ndarray _spupil   # small pupil (without guard band)
    cdef np.ndarray _apodizer # apodizer (same size as small pupil)
    cdef long  _p1         # min x,y for valid points in mpupil
    cdef long  _p2         # max x,y for valid points in mpupil
    cdef long  _n          # linear size of mpupil
    cdef long  _n1         # min x,y for valid points in ipupil
    cdef long  _n2         # max x,y for valid points in ipupil


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
    cdef float zerop
    """target flux for magnitude 0"""
    cdef np.ndarray  dms_seen
    """index of dms seen by the target"""


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
    cdef readonly np.ndarray dms_seen
    """index of dms seen by the WFS"""


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
    cdef readonly long  _nvalid         # number of valid subaps

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

    cdef readonly np.ndarray _lgskern      # lgs kernels for each subap
    cdef readonly np.ndarray _profna       # sodium profile
    cdef readonly np.ndarray _altna        # corresponding altitude
    cdef readonly np.ndarray _prof1d       # hr profile
    cdef readonly np.ndarray _profcum      # hr profile cumulated
    cdef readonly np.ndarray _beam         # 1d beam function
    cdef readonly np.ndarray _ftbeam       # 1d beam function fft
    cdef readonly np.ndarray _azimuth      # angles of rotation for each spot

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



    cpdef prep_lgs_prof(self,int nsensors, Param_tel p_tel, 
                        np.ndarray[dtype=np.float32_t] prof,
                        np.ndarray[dtype=np.float32_t] h,
                        float beam, Sensors sensors,
                        bytes center=?, int imat=?)

    cdef make_lgs_prof1d(self, Param_tel p_tel,
            np.ndarray[dtype=np.float32_t] prof, 
            np.ndarray[dtype=np.float32_t] h,
            float beam, bytes center=?)

#################################################
# P-Class (parametres) Param_dm
#################################################
cdef class Param_dm:
  cdef readonly bytes  type_dm
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
  cdef readonly np.ndarray _influstart
  cdef readonly np.ndarray _influkernel # Array to convolution kernel for comp_dmshape
  cdef readonly list _klbas
  """ np.ndarray to a kl struct"""


  #cdef set_xpos(self,np.ndarray[ndim=1,dtype=np.float32_t] xpos)
  #cdef set_ypos(self,np.ndarray[ndim=1,dtype=np.float32_t] ypos)

#################################################
# P-Class (parametres) Param_rtc
#################################################
cdef class Param_rtc:
    cdef long    nwfs           # number of wfs
    cdef list centroiders     # an array of centroiders
    cdef list controllers     #an array of controllers


#################################################
# P-Class (parametres) Param_centroider
#################################################
cdef class Param_centroider:
    cdef long    nwfs       # index of wfs in y_wfs structure on which we want to do centroiding
    cdef bytes  type_centro # type of centroiding "cog", "tcog", "bpcog", "wcog", "corr"
    cdef bytes  type_fct    # type of ref function "gauss", "file", "model"
    cdef np.ndarray weights    # optional reference function(s) used for centroiding
    cdef long    nmax       #
    cdef float   thresh     #
    cdef float   width      # width of the Gaussian
    cdef long    sizex      #
    cdef long    sizey      #
    cdef np.ndarray interpmat  # optional reference function(s) used for centroiding


#################################################
# P-Class (parametres) Param_controller
#################################################
cdef class Param_controller:
    cdef readonly bytes  type_control   # type of controller
    cdef readonly np.ndarray nwfs        # index of wfss in controller
    cdef readonly np.ndarray nvalid      # number of valid subaps per wfs
    cdef readonly np.ndarray ndm         # index of dms in controller
    cdef readonly np.ndarray nactu       # number of controled actuator per dm
    cdef readonly np.ndarray imat        # full interaction matrix
    cdef readonly np.ndarray cmat        # full control matrix
    cdef readonly float   maxcond        # max condition number
    cdef float    delay         # 
    cdef float   gain           #
    cdef long    nkl            # Florain features : number of KL modes used for computation of covmat in case of minimum variance controller
    cdef long    cured_ndivs    # subdivision levels in cured
    cdef int     modopti        # Flag for modal optimization
    cdef int     nrec           # Number of sample of open loop slopes for modal optimization computation
    cdef int     nmodes         # Number of modes for M2V matrix (modal optimization)
    cdef float     gmin         # Minimum gain for modal optimization
    cdef float     gmax         # Maximum gain for modal optimization
    cdef int     ngain          # Number of tested gains

#################################################
# P-Class (parametres) Param_loop
#################################################
cdef class Param_loop:
    cdef readonly long niter     # number of iterations
    cdef readonly float ittime   # iteration time (in sec)


#################################################
# P-Class (parametres) Param_kl_basis
#################################################
cdef class Param_kl_basis:
    cdef long _nr
    cdef long _ni
    cdef long _np
    cdef long _nfunct
    cdef float _cobs
    cdef long _nord
    cdef np.ndarray _radp
    cdef np.ndarray _evals
    cdef np.ndarray _npo
    cdef np.ndarray _ord
    cdef np.ndarray _rabas
    cdef np.ndarray _azbas
    #######################
    cdef long _ncp
    cdef long _ncmar
    cdef np.ndarray _px
    cdef np.ndarray _py
    cdef np.ndarray _cr
    cdef np.ndarray _cp
    cdef np.ndarray _pincx
    cdef np.ndarray _pincy
    cdef np.ndarray _pincw
    cdef np.ndarray _ap
