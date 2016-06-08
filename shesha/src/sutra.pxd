include "../par.pxi"

from naga_context cimport carma_context, carma_device
from naga_obj cimport carma_obj
from naga_host_obj cimport carma_host_obj
from naga_sparse_obj cimport carma_sparse_obj
from naga_streams cimport carma_streams
from naga_obj cimport cuFloatComplex, cuDoubleComplex, cufftHandle, tuple_t
# from naga cimport *

cimport numpy as np

from libc.stdlib cimport malloc, free
from libcpp.cast cimport dynamic_cast
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool

from libc.stdint cimport uintptr_t

from cpython.string cimport PyString_AsString
from libc.math cimport sin

cdef np.float32_t dtor = (np.pi / 180.)


ctypedef enum cudaChannelFormatDesc:
    f
    w
    x
    y
    z


cdef extern from "sutra_ao_utils.h":
    int roll_mult[float](float * odata, float * idata, int N, int M, float alpha, carma_device * device)

#################################################
# C-Class sutra_phase
#################################################
cdef extern from "sutra_phase.h":
    cdef cppclass sutra_phase:
        carma_obj[float] * d_screen
        long screen_size
        float * zerCoeff
        carma_obj[float] * zernikes
        carma_obj[float] * mat

#################################################
# C-Class sutra_telescope
#################################################
cdef extern from "sutra_telescope.h":
    cdef cppclass sutra_telescope:

        carma_context * current_context
        int device  # // device

        long pup_size  # // size of pupil
        long num_eleme_pup  # // number of points in the pupil

        carma_obj[float] * d_pupil  # // the pupil mask
        carma_obj[float] * d_phase_ab_M1  # // the phase aberration for M1

        long pup_size_m  # // size of pupil

        carma_obj[float] * d_pupil_m  # // the pupil mask
        carma_obj[float] * d_phase_ab_M1_m  # // the phase aberration for M1

        sutra_telescope(carma_context * context, long pup_size, long num_eleme_pup, float * pupil, float * phase_ab_M1, long pup_size_m, float * pupil_m, float * phase_ab_m1_m)


cdef extern from "sutra_turbu.h":
    #################################################
    # C-Class sutra_t_screen
    #################################################
    cdef cppclass sutra_tscreen:
        int device  # The device #
        sutra_phase * d_tscreen  # The phase screen
        long screen_size  # size of phase screens
        float amplitude  # amplitude for extrusion (r0^-5/6)
        float altitude
        float windspeed
        float winddir
        float deltax  # number of rows to extrude per iteration
        float deltay  # number of lines to extrude per iteration
        # internal
        float accumx
        float accumy

        float norm_vk
        bool vk_on
        carma_context * current_context

        carma_obj[float] * d_A  # A matrix for extrusion
        carma_obj[float] * d_B  # B matrix for extrusion
        carma_obj[unsigned int] * d_istencilx  # stencil for column extrusion
        carma_obj[unsigned int] * d_istencily  # stencil for row extrusion
        carma_obj[float] * d_z  # tmp array for extrusion process
        carma_obj[float] * d_noise  # tmp array containing random numbers
        carma_obj[float] * d_ytmp

        sutra_tscreen(carma_context * context, long size, long size2, float amplitude,
                      float altitude, float windspeed, float winddir, float deltax,
                      float deltay, int device)
        # sutra_tscreen(const sutra_tscreen& tscreen)

        int init_screen(float * h_A, float * h_B, unsigned int * h_istencilx,
                        unsigned int * h_istencily, int seed)
        int extrude(int dir)
        int init_vk(int seed, int pupd)
        int generate_vk(float l0, int nalias)


#################################################
# C-Class sutra_atmos
#################################################
    cdef cppclass sutra_atmos:
        int nscreens
        map.map[float, sutra_tscreen * ] d_screens
        float r0
        # carma_obj[float] *d_pupil
        carma_context * current_context

        sutra_atmos(carma_context * context, int nscreens, float * r0, long * size,
                    long * size2, float * altitude, float * 
                    windspeed, float * winddir,
                    float * deltax, float * deltay, int device)
        init_screen(float alt, float * h_A, float * h_B, unsigned int * h_istencilx,
                    unsigned int * h_istencily, int seed)
        int move_atmos()


cdef extern from "sutra_target.h":
    #################################################
    # C-Class sutra_source
    #################################################
    cdef cppclass sutra_source:

        int device  # device #
        float tposx  # x position of target on the sky
        float tposy  # y position of target on the sky
        long npos  # number of points in the pupil
        float mag  # brightness of target
        float Lambda "lambda"  # imaging lambda
        float zp  # imaging zero point
        float scale  # phase scale
        bool lgs  # flag for lgs
        char * type  # type of source : target / wfs
        int block_size  # optimum block size of device

        float strehl_se  # short exposure strehl
        float strehl_le  # long exposure strehl
        float ref_strehl  # reference for strehl computation
        int strehl_counter  # counter for le strehl computation
        float phase_var  # current phase variance in the pupil
        float phase_var_avg  # average phase variance in the pupil
        int phase_var_count  # counter for average phase variance in the pupil

        sutra_phase * d_phase  # phase for this target
        # INTRO PHASE INSTRU
        # sutra_phase *d_phase_instru;
        #
        carma_host_obj[float] * phase_telemetry  #
        sutra_lgs * d_lgs  # the lgs object
        carma_obj[float] * object  # the object intensity map
        carma_obj[float] * d_pupil  # the pupil mask
        carma_obj[float] * d_image  # the resulting image for this target
        carma_obj[float] * d_leimage  # the long exposure image for this target
        # the complex amplitude in the pupil plane
        carma_obj[cuFloatComplex] * d_amplipup
        # the valid phase points in the pupil (target only)
        carma_obj[float] * d_phasepts
        # positions of valid phase points in the pupil (target only)
        carma_obj[int] * d_wherephase
        # type_screen unknown cf sutra_DM
        map[pair[string, int], float] xoff  # x reference for raytracing
        map[pair[string, int], float] yoff  # y reference for raytracing
        carma_context * current_context

        int add_layer(char * l_type, float alt, float xoff, float yoff)
        int init_strehlmeter()
        int raytrace(sutra_atmos * atmos)
        int raytrace(sutra_dms * ydms, int rst, int do_phase_var)
        int comp_image(int puponly, bool comp_le)
        int comp_strehl()

#################################################
# C-Class sutra_target
#################################################
    cdef cppclass sutra_target:
        int ntargets
        vector.vector[sutra_source * ] d_targets

        sutra_target(carma_context * context, sutra_telescope * d_tel, int ntargets, float * xpos, float * ypos,
                     float * Lambda, float * mag, float zerop, long * sizes, int Npts, int device)

        # not implemented
        # int get_phase(int ntarget, float *dest)

    cdef int target_texraytrace(float * d_odata, float * d_idata, int nx, int ny, int Nx,
                                int Ny, float xoff, float yoff, int Ntot, cudaChannelFormatDesc channelDesc,
                                carma_device * device)

    cdef int target_raytrace(float * d_odata, float * d_idata, int nx, int ny, int Nx,
                             float xoff, float yoff, int block_size)

    cdef int target_lgs_raytrace(float * d_odata, float * d_idata, int nx, int ny, int Nx,
                                 float xoff, float yoff, float delta, int block_size)

    cdef int target_raytrace_async(carma_streams streams, float * d_odata, float * d_idata,
                                   int nx, int ny, int Nx, float xoff, float yoff, int block_size)

    cdef int target_raytrace_async(carma_host_obj[float] * phase_telemetry, float * d_odata,
                                   float * d_idata, int nx, int ny, int Nx, float xoff, float yoff, int block_size)

    # cdef int fft_goodsize(long size)


cdef extern from "sutra_wfs.h":
    #################################################
    # C-Class sutra_sensor
    #################################################
    cdef cppclass sutra_sensors:
        int device
        bool error_budget
        carma_context * current_context
        size_t nsensors()

        vector[sutra_wfs * ] d_wfs
        map[vector[int], cufftHandle * ] campli_plans
        map[vector[int], cufftHandle * ] fttotim_plans
        map[vector[int], cufftHandle * ] ftlgskern_plans

        carma_obj[cuFloatComplex] * d_camplipup
        carma_obj[cuFloatComplex] * d_camplifoc
        carma_obj[cuFloatComplex] * d_fttotim
        carma_obj[cuFloatComplex] * d_ftlgskern
        carma_obj[float] * d_lgskern

        sutra_sensors(carma_context * context, sutra_telescope * d_tel, const char ** type, int nwfs, long * nxsub,
                      long * nvalid, long * npix, long * 
                      nphase, long * nrebin, long * nfft,
                      long * ntot, long * npup, float * pdiam, float * nphot, float * nphot4imat, int * lgs, int device, bool error_budget)
        sutra_sensors(carma_context * context, sutra_telescope * d_tel, int nwfs, long * nxsub, long * nvalid,
                      long * nphase, long npup, float * pdiam, int device)

        int sensors_initgs(float * xpos, float * ypos, float * Lambda, float * mag, float zerop,
                           long * size, float * noise, long * seed)
        int sensors_initgs(float * xpos, float * ypos, float * Lambda, float * mag, float zerop,
                           long * size, float * noise)
        int sensors_initgs(float * xpos, float * ypos, float * Lambda, float * mag, float zerop,
                           long * size)
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
        float nphot4imat
        float noise
        bool lgs
        bool kernconv
        bool error_budget

        int rank
        int offset
        int nvalid_tot
        int * displ_bincube
        int * count_bincube

        cufftHandle * campli_plan
        cufftHandle * fttotim_plan
        carma_obj[cuFloatComplex] * d_ftkernel
        carma_obj[cuFloatComplex] * d_camplipup
        carma_obj[cuFloatComplex] * d_camplifoc
        carma_obj[cuFloatComplex] * d_fttotim

        carma_obj[float] * d_pupil
        carma_obj[float] * d_bincube
        carma_obj[float] * d_bincube_notnoisy
        carma_obj[float] * d_binimg
        carma_obj[float] * d_subsum
        carma_obj[float] * d_offsets
        carma_obj[float] * d_fluxPerSub
        carma_obj[float] * d_sincar
        carma_obj[int] * d_hrmap

        carma_obj[float] * d_slopes

        carma_host_obj[float] * image_telemetry

      # sh only
        carma_obj[int] * d_phasemap
        carma_obj[int] * d_binmap
        carma_obj[int] * d_validsubsx  # nvalid
        carma_obj[int] * d_validsubsy  # nvalid
        carma_obj[int] * d_istart  # nxsub
        carma_obj[int] * d_jstart  # nxsub

      # pyramid only
        carma_obj[float] * d_hrimg
        carma_obj[float] * d_submask
        carma_obj[float] * d_psum
        carma_obj[cuFloatComplex] * d_phalfxy
        carma_obj[cuFloatComplex] * d_poffsets

        carma_host_obj[float] * pyr_cx
        carma_host_obj[float] * pyr_cy

        sutra_source * d_gs

        carma_streams * streams
        int nstreams

        carma_context * current_context

        sutra_wfs(carma_context * context, sutra_telescope * d_tel, sutra_sensors * sensors, const char * type, long nxsub, long nvalid,
                  long npix, long nphase, long nrebin, long nfft, long ntot, long npup,
                  float pdiam, float nphotons, float nphot4imat, int lgs, int device)
        sutra_wfs(carma_context * context, sutra_telescope * d_tel, long nxsub, long nvalid, long nphase,
                  long npup, float pdiam, int device)
        sutra_wfs(const sutra_wfs & wfs)

        int wfs_initarrays(int * phasemap, int * hrmap, int * binmap, float * offsets,
                           float * fluxPerSub, int * validsubsx,
                           int * validsubsy, int * istart, int * jstart, cuFloatComplex * kernel)
        int wfs_initarrays(int * phasemap, float * offsets, float * fluxPerSub,
                           int * validsubsx, int * validsubsy, int * istart, int * jstart)
        int wfs_initarrays(cuFloatComplex * halfxy, cuFloatComplex * offsets,
                           float * focmask, float * cx, float * cy,
                           float * sincar, int * phasemap, int * validsubsx, int * validsubsy)
        int wfs_initgs(sutra_sensors * sensors, float xpos, float ypos, float Lambda, float mag, long size,
                       float noise, long seed)
        int load_kernels(float * lgskern)
        int sensor_trace(sutra_atmos * yatmos)
        int sensor_trace(sutra_atmos * atmos, sutra_dms * ydms)
        int sensor_trace(sutra_dms * ydm, int rst)
        int comp_image()
        int comp_generic()
        int define_mpi_rank(int rank, int size)
        int allocate_buffers(sutra_sensors * sensors)
        int fill_binimage(int async)


#################################################
# C-Class sutra_wfs_sh
#################################################
cdef extern from "sutra_wfs_sh.h":
    cdef cppclass sutra_wfs_sh(sutra_wfs):
        int  wfs_initarrays(int * phasemap, int * hrmap, int * binmap, float * offsets,
                            float * fluxPerSub, int * validsubsx,
                            int * validsubsy, int * istart, int * jstart, cuFloatComplex * kernel)
        int slopes_geom(int type, float * slopes)
        int slopes_geom(int type)


#################################################
# C-Class sutra_wfs_pyr_roof
#################################################
cdef extern from "sutra_wfs_pyr.h":
    cdef cppclass sutra_wfs_pyr(sutra_wfs):
        int wfs_initarrays(cuFloatComplex * halfxy, cuFloatComplex * offsets,
                           float * focmask, float * cx, float * cy,
                           float * sincar, int * phasemap, int * validsubsx, int * validsubsy)

#################################################
# C-Class sutra_wfs_pyr_roof
#################################################
cdef extern from "sutra_wfs_pyr_roof.h":
    cdef cppclass sutra_wfs_pyr_roof(sutra_wfs_pyr):
        sutra_wfs_pyr_roof(const sutra_wfs_pyr_roof & wfs)

#################################################
# C-Class sutra_wfs_pyr_pyr4
#################################################
cdef extern from "sutra_wfs_pyr_pyr4.h":
    cdef cppclass sutra_wfs_pyr_pyr4(sutra_wfs_pyr):
        sutra_wfs_pyr_pyr4(const sutra_wfs_pyr_pyr4 & wfs)

#################################################
# C-Class sutra_wfs_pyr_pyrhr
#################################################
cdef extern from "sutra_wfs_pyr_pyrhr.h":
    cdef cppclass sutra_wfs_pyr_pyrhr(sutra_wfs_pyr):
        sutra_wfs_pyr_pyrhr(const sutra_wfs_pyr_pyrhr & wfs)
        int wfs_initarrays(cuFloatComplex * halfxy, float * cx, float * cy,
                           float * sincar, int * validsubsx, int * validsubsy, int *phasemap, float *fluxPerSub)
        void comp_modulation(int cpt);
        int slopes_geom(int type, float * slopes)
        int slopes_geom(int type)

#################################################
# C-Class sutra_wfs_geom
#################################################
cdef extern from "sutra_wfs_geom.h":
    cdef cppclass sutra_wfs_geom(sutra_wfs):
        # sutra_wfs_geom(carma_context *context, long nxsub, long nvalid, long nphase,
        #    long npup, float pdiam, int device)
        sutra_wfs_geom(const sutra_wfs_geom & wfs)
        int wfs_initarrays(int * phasemap, float * offsets,
                           float * fluxPerSub, int * validsubsx, int * validsubsy)


#################################################
# C-Class sutra_centroider
#################################################
cdef extern from "sutra_centroider.h":
    cdef cppclass sutra_centroider:
        int device
        sutra_wfs * wfs
        int nwfs
        int nvalid

        float offset
        float scale

        carma_context * current_context

        bool is_type(string typec)

        string get_type()

        int get_cog(carma_streams * streams, float * cube, float * subsum, float * centroids,
                    int nvalid, int npix, int ntot)
        int get_cog(float * subsum, float * slopes)
        int get_cog()

    int convert_centro(float * d_odata, float * d_idata, float offset, float scale,
                       int N, carma_device * device)


#################################################
# C-Class sutra_centroider_tcog
#################################################
cdef extern from "sutra_centroider_tcog.h":
    cdef cppclass sutra_centroider_tcog(sutra_centroider):
        int set_threshold(float threshold)
        bool is_type(string typec)
        # string get_type()

#################################################
# C-Class sutra_centroider_corr
#################################################
cdef extern from "sutra_centroider_corr.h":
    cdef cppclass sutra_centroider_corr(sutra_centroider):
        int init_bincube()
        bool is_type(string typec)
        string get_type()

        int init_corr(int isizex, int isizey, float * interpmat)
        int load_corr(float * corr, float * corr_norm, int ndim)


#################################################
# C-Class sutra_centroider_wcog
#################################################
cdef extern from "sutra_centroider_wcog.h":
    cdef cppclass sutra_centroider_wcog(sutra_centroider):
        int npix
        carma_obj[float] * d_weights

        sutra_centroider_wcog(carma_context * context, sutra_sensors * sensors,
                              int nwfs, long nvalid, float offset, float scale, int device)
        sutra_centroider_wcog(const sutra_centroider & centroider)
        string get_type()
        int init_weights()
        int load_weights(float * weights, int ndim)
        int get_cog(carma_streams * streams, float * cube, float * subsum,
                    float * centroids, int nvalid, int npix, int ntot)
        int get_cog(float * subsum, float * slopes)
        int get_cog()


'''
#################################################
# Dynamic casts
#################################################
cdef extern from *:
    sutra_wfs_geom* dynamic_cast_wfs_geom_ptr "dynamic_cast<sutra_wfs_geom*>" (sutra_wfs*) except NULL
    sutra_wfs_sh* dynamic_cast_wfs_sh_ptr "dynamic_cast<sutra_wfs_sh*>" (sutra_wfs*) except NULL
    sutra_wfs_pyr_pyr4* dynamic_cast_wfs_pyr_pyr4_ptr "dynamic_cast<sutra_wfs_pyr_pyr4*>" (sutra_wfs*) except NULL
    sutra_wfs_pyr_pyrhr* dynamic_cast_wfs_pyr_pyrhr_ptr "dynamic_cast<sutra_wfs_pyr_pyrhr*>" (sutra_wfs*) except NULL
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
'''

cdef extern from "sutra_dm.h":
    #################################################
    # C-Class sutra_dms
    #################################################
    ctypedef pair[string, float] type_screen

    cdef cppclass sutra_dms:
        vector[sutra_dm * ] d_dms
        vector[type_screen] d_type

        sutra_dms(int ndm)

        int add_dm(carma_context * context, const char * type, float alt, long dim,
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

        sutra_phase * d_shape

        carma_obj[float] * d_comm

        carma_obj[float] * d_influ  # if relevant

        # pzt
        carma_obj[int] * d_influpos
        carma_obj[int] * d_npoints
        carma_obj[int] * d_istart
        carma_obj[int] * d_xoff
        carma_obj[int] * d_yoff
        carma_obj[int] * d_pos  # Convolution preprocess
        carma_obj[float] * d_KLbasis
        carma_obj[float] * d_kernconv  # Convolution preprocess
        carma_obj[float] * d_mapactu  # Convolution process
        carma_obj[cuFloatComplex] * d_ftkernconv  # Convolution process
        carma_obj[cuFloatComplex] * d_ftmapactu  # Convolution process

        # zernike
        carma_obj[float] * d_coeffs
        carma_obj[float] * d_mask
        carma_obj[float] * d_zr
        carma_obj[float] * d_ztheta

        # sutra_kl *d_kl

        carma_context * current_context

        int nact()
        int pzt_loadarrays(float * influ, float * influ2, tuple_t[float] * influ3,
                           int * influpos, int * influpos2, int * 
                           npoints, int * istart,
                           int * xoff, int * yoff, float * kernconv)
        int kl_loadarrays(float * rabas, float * azbas, int * ord, float * cr, float * cp)

        int reset_shape()

        int comp_oneactu(int nactu, float ampli)

        int comp_shape()
        int comp_shape(float * comm)
        int compute_KLbasis(float * xpos, float * ypos, int * indx, long dim,
                            float norm, float ampli)

#################################################
# C-Class sutra_controller
#################################################
cdef extern from "sutra_controller.h":

    cdef cppclass sutra_controller:
        # allocation of d_centroids and d_com
        sutra_controller(carma_context * context, int nslope, int nactu, float delay,
                         sutra_dms * dms, char ** type, float * alt, int ndm)

        string get_type()

#!!!! YOU MUST set d_centroids before calling it!!!!
        int comp_com()

# It is better to have something like this (+protected d_centroids):
# virtual int comp_com (carma_obj<float> *new_centroids)=0
# it would imply copy, but would be much safer

        inline int nactu()
        inline int nslope()

        # TODO cublasHandle_t cublas_handle()

        int set_perturbcom(float * perturb, int N)
        int set_openloop(int open_loop_status)
        int comp_voltage()
        int syevd_f(char meth, carma_obj[float] * d_U, carma_host_obj[float] * h_eingenvals)
        int invgen(carma_obj[float] * d_mat, float cond, int job)
        int command_delay()

        # I would propose to make them protected (+ proper
        # set of fuctions). It could make life easier!
        # But we should discuss it
        int cpt_pertu
        int open_loop
        float delay
        # Coefficient for linear interpolation on command buffer to allow
        # non-integer delay
        float a
        # Coefficient for linear interpolation on command buffer to allow
        # non-integer delay
        float b
        # Coefficient for linear interpolation on command buffer to allow
        # non-integer delay
        float c
        vector[sutra_dm * ] d_dmseen
        carma_obj[float] * d_centroids  # current centroids
        carma_obj[float] * d_com  # current command
        carma_obj[float] * d_perturb  # perturbation command buffer
        carma_obj[float] * d_voltage  # commands sent to mirror
        carma_obj[float] * d_com1  # commands k-1
        carma_obj[float] * d_com2  # commands k-2

        carma_streams * streams

        int device
        carma_context * current_context


#################################################
# C-Class sutra_controller_generic
#################################################
cdef extern from "sutra_controller_generic.h":
    cdef cppclass sutra_controller_generic(sutra_controller):
        carma_obj[float] * d_matE
        carma_obj[float] * d_cmat
        carma_obj[float] * d_gain
        carma_obj[float] * d_decayFactor
        carma_obj[float] * d_compbuff
        string command_law

        sutra_controller_generic(carma_context * context, long nvalid, long nactu,
                                 float delay, sutra_dms * dms, char ** type, float * alt, int ndm)
        sutra_controller_generic(const sutra_controller_generic & controller)

        string get_type()
        string get_commandlaw()
        int set_decayFactor(float * decayFactor)
        int set_mgain(float * gain)
        int set_cmat(float * cmat)
        int set_matE(float * matE)
        int set_commandlaw(string law)
        int comp_com()


#################################################
# C-Class sutra_controller_geo
#################################################
cdef extern from "sutra_controller_geo.h":
    cdef cppclass sutra_controller_geo(sutra_controller):
        float gain
        long Nphi

        carma_obj[float] * d_gain
        carma_obj[float] * d_proj
        carma_obj[double] * d_phi
        carma_obj[int] * d_indx_pup
        carma_sparse_obj[double] * d_IFsparse
        carma_obj[float] * d_geocov
        carma_obj[double] * d_compdouble
        carma_obj[float] * d_compfloat

        sutra_controller_geo(carma_context * context, long nactu, long Nphi,
                             float delay, sutra_dms * dms, char ** type, float * alt, int ndm, bool wfs_direction)
        sutra_controller_geo(const sutra_controller_geo & controller)

        string get_type()

        int load_Btt(float *Btt)
        int set_gain(float gain)
        int load_mgain(float * mgain)
        int comp_dphi(sutra_source * target, bool wfs_direction)
        int comp_com()
        int init_proj(sutra_dms * dms, int * indx_dm, float * unitpervolt, int * indx_pup)
        int init_proj_sparse(sutra_dms * dms, int * indx_dm, float * unitpervolt, int * indx_pup, int * indx_mpup)


#################################################
# C-Class sutra_controller_ls
#################################################
cdef extern from "sutra_controller_ls.h":
    cdef cppclass sutra_controller_ls(sutra_controller):
        float gain

        carma_obj[float] * d_imat
        carma_obj[float] * d_cmat
        carma_obj[float] * d_gain

# svd computations
        carma_obj[float] * d_eigenvals
        carma_host_obj[float] * h_eigenvals
        carma_obj[float] * d_U

# loop components
        carma_obj[float] * d_cenbuff  # centroids circular buffer
        carma_obj[float] * d_err  # current error

# Modal optimization components
        int is_modopti  # Flag for using modal optimization
        int nrec  # Number of recorded open slopes measurements
        int nmodes  # Number of modes
        float gmin  # Gain min
        float gmax  # Gain max
        int ngain  # Number of tested gains between gmin and gmax
        float Fs  # Sampling frequency
        int cpt_rec  # Counter for modal gains refresh
        carma_obj[float] * d_M2V  # Modes to Volts matrix
        carma_obj[float] * d_S2M  # Slopes to Modes matrix
        # Open-loop measurements buffer, recorded and loaded from Yorick
        carma_obj[float] * d_slpol
        carma_obj[float] * d_Hcor  # Transfer function
        carma_obj[float] * d_com1  # Command k-1 for POLC
        carma_obj[float] * d_com2  # Command k-2 for POLC
        carma_obj[float] * d_compbuff  # Buffer for POLC computation
        carma_obj[float] * d_compbuff2  # Buffer for POLC computation

        sutra_controller_ls(carma_context * context, long nvalid, long nactu,
                            float delay, sutra_dms * dms, char ** type, float * alt, int ndm)
        sutra_controller_ls(const sutra_controller_ls & controller)

        string get_type()

        int svdec_imat()
        int build_cmat(int nfilt, bool filt_tt)
        int build_cmat(int nfilt)
        int build_cmat_modopti()
        int frame_delay()
        int comp_com()
        int set_gain(float gain)
        int set_mgain(float * mgain)
        int set_cmat(float * cmat)
        int set_delay(float delay)
        int init_modalOpti(int nmodes, int nrec, float * M2V, float gmin, float gmax, int ngain, float Fs)
        int loadOpenLoopSlp(float * ol_slopes)
        int modalControlOptimization()
        int compute_Hcor()


#################################################
# C-Class sutra_controller_mv
#################################################
cdef extern from "sutra_controller_mv.h":
    cdef cppclass sutra_controller_mv(sutra_controller):

        float gain

        carma_obj[float] * d_imat
        carma_obj[float] * d_cmat
        carma_obj[float] * d_gain

# Cphim & Cmm features
        carma_obj[float] * d_covmat
        carma_obj[float] * d_KLbasis
        carma_obj[float] * d_noisemat
        carma_obj[float] * d_Cmm
        carma_obj[float] * d_Cphim
# svd computations
        carma_host_obj[float] * h_Cmmeigenvals
        carma_host_obj[float] * h_eigenvals
# carma_obj[float] *d_U

# loop components
        carma_obj[float] * d_cenbuff  # centroids circular buffer
        carma_obj[float] * d_com1  # commands k-1 (for POLC)
        carma_obj[float] * d_com2  # commands k-2 (for POLC)
        carma_obj[float] * d_compbuff  # Buffer for computations
        carma_obj[float] * d_compbuff2
        carma_obj[float] * d_olmeas  # Open-loop measurements for POLC
        carma_obj[float] * d_err  # current error

        # TODO cublasHandle_t cublas_handle

        sutra_controller_mv(carma_context * context, long nvalid, long nactu,
                            float delay, sutra_dms * dms, char ** type, float * alt, int ndm)
        sutra_controller_mv(const sutra_controller_mv & controller)

        string get_type()

        int svdec_imat()
        int build_cmat(const char * dmtype, char * method)
        int build_cmat(float cond)
        int frame_delay()
        int comp_com()
        int set_gain(float gain)
        int set_mgain(float * mgain)
        int set_delay(float delay)
        int set_cmat(float * cmat)
        # Florian features
        int load_noisemat(float * noise)
        int do_covmat(sutra_dm * ydm, char * method, int * indx_pup, long dim, float * xpos, float * ypos, long Nkl, float norm, float ampli)
        int do_geomat(carma_obj[float] * d_geocov, carma_obj[float] * d_IF, long n_pts, float ampli)
        int piston_filt(carma_obj[float] * d_statcov)
        int piston_filt_cphim(carma_obj[float] * d_cphim, float * F)
        int filter_cphim(float * F, float * Nact)
        int filter_cmat(float cond)
        int invgen(carma_obj[float] * d_mat, float cond, int job)
        int invgen(carma_obj[float] * d_mat, carma_host_obj[float] * h_eigen, float cond)
        int invgen_cpu(carma_obj[float] * d_mat, carma_host_obj[float] * h_eigen, float cond)
        # int do_statmat(float *statcov,long dim, float *xpos, float *ypos,
        # float norm, carma_device *device)
        int DDiago(carma_obj[float] * d_statcov, carma_obj[float] * d_geocov)
        int load_covmat(float * covmat)
        int load_klbasis(float * klbasis)
        int compute_Cmm(sutra_atmos * atmos, sutra_sensors * sensors, double * L0, double * cn2, double * alphaX, double * alphaY, double diamTel, double cobs)
        int compute_Cphim(sutra_atmos * atmos,
                          sutra_sensors * sensors, sutra_dms * 
                          dms, double * L0, double * cn2,
                          double * alphaX, double * alphaY, double * 
                          X, double * Y, double * xactu,
                          double * yactu, double diamTel, double * k2, long * NlayerDm,
                          long * indLayerDm, double FoV, double * pitch, double * alt_dm)


#################################################
# C-Class sutra_controller_cured
#################################################
cdef extern from "sutra_controller_cured.h":
    cdef cppclass sutra_controller_cured(sutra_controller):
        float gain
        int   ndivs  # number of subdivision levels for cured
        bool  tt_flag  # flag for separate tt

        # data for CuReD */
        carma_host_obj[float] * h_centroids
        carma_host_obj[float] * h_err
        carma_obj[float] * d_err  # current error
        carma_obj[float] * d_cenbuff  # centroids circular buffer

        # data for CuReD */
        carma_obj[float] * d_imat

        # structures needed to run CuReD */
        # sysCure* h_syscure
        void * h_syscure
        # parCure* h_parcure
        void * h_parcure

        sutra_controller_cured(carma_context * context, long nvalid, long nactu,
                               float delay, sutra_dms * dms, char ** type, float * alt, int ndm)
        sutra_controller_cured(const sutra_controller_cured & controller)

        string get_type()
        int set_gain(float gain)

        int comp_com()

        int init_cured(int nxsubs, int * isvalid, int ndivs, int tt)
        int frame_delay()
        int set_delay(float delay)


#################################################
# C-Class sutra_controller_kalman
#################################################
cdef extern from "sutra_controller_kalman.h":
    cdef cppclass sutra_controller_kalman(sutra_controller):
        sutra_controller_kalman(carma_context * context, int nvalid_, int nactu_, sutra_dms * dms, char ** type, float * alt, int ndm)

        void init_kalman(carma_host_obj[float] & D_Mo, carma_host_obj[float] & N_Act, carma_host_obj[float] & PROJ, bool is_zonal, bool is_sparse, bool is_GPU)

        double gettime()
        double gettime_op1()
        double gettime_op2()
        double gettime_op3()

        void calculate_gain(float bruit, carma_host_obj[float] & SigmaV,
                            carma_host_obj[float] & atur, carma_host_obj[float] & btur)

        int set_gain(float k_W)

        # TODO cusparseHandle_t cusparseHandle
        # kp_kalman_core_sparse* core_sparse
        # kp_kalman_core_full* core_full
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

        vector[sutra_centroider * ] d_centro
        vector[sutra_controller * ] d_control

        carma_context * current_context

        sutra_rtc(carma_context * context)
        # sutra_rtc(const sutra_rtc& yrtc)

        int add_centroider(sutra_sensors * sensors, int nwfs, long nvalid, float offset, float scale, long device,
                           char * typec)
        int rm_centroider()
        int add_controller_geo(int nactu, int Nphi, float delay, long device,
                               sutra_dms * dms, char ** type_dmseen, float * alt, int ndm, bool wfs_direction)
        int add_controller(int nactu, float delay, long device, const char * typec,
                           sutra_dms * dms, char ** type_dmseen, float * alt, int ndm)
        int rm_controller()

        int do_imat(int ncntrl, sutra_dms * ydms)
        # int do_imatkl(int ncntrl, sutra_dms *ydms)
        # int do_imatkl4pzt(int ncntrl, sutra_dms *ydms)
        int do_imat_geom(int ncntrl, sutra_dms * ydm, int type)

        int do_centroids()
        int do_centroids(int ncntrl)
        int do_centroids(int ncntrl, bool imat)
        int do_centroids_geom(int ncntrl)
        int do_control(int ncntrl)
        int apply_control(int ncntrl, sutra_dms * ydm)

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

        cufftHandle * ftlgskern_plan
        carma_obj[float] * d_doffaxis
        carma_obj[float] * d_azimuth
        carma_obj[float] * d_prof1d
        carma_obj[float] * d_profcum
        carma_obj[cuFloatComplex] * d_prof2d
        carma_obj[float] * d_beam
        carma_obj[cuFloatComplex] * d_ftbeam
        carma_obj[float] * d_lgskern
        carma_obj[cuFloatComplex] * d_ftlgskern

        carma_context * current_context

        sutra_lgs(carma_context * context, sutra_sensors * sensors, long nvalid, long npix, long nmaxhr)
        sutra_lgs(const sutra_lgs & lgs)

        int lgs_init(int nprof, float hg, float h0, float deltah, float pixsie,
                     float * doffaxis, float * 
                     prof1d, float * profcum, float * beam,
                     cuFloatComplex * ftbeam, float * azimuth)
        int load_prof(float * prof1d, float * profcum, float hg, float h0, float deltah)
        int lgs_update(carma_device * device)
        int lgs_makespot(carma_device * device, int nin)
        int load_kernels(float * lgskern, carma_device * device)

IF USE_BRAMA == 1:
    #################################################
    # C-Class sutra_rtc_brama
    #################################################
    cdef extern from "sutra_rtc_brama.h":
        cdef cppclass sutra_rtc_brama(sutra_rtc):

            sutra_rtc_brama(
                carma_context * context, sutra_sensors * wfs, sutra_target * target, char * name)

            void publish()

    #################################################
    # C-Class sutra_target_brama
    #################################################
    cdef extern from "sutra_target_brama.h":
        cdef cppclass sutra_target_brama(sutra_target):

            sutra_target_brama(carma_context * context, char * name, sutra_telescope * d_tel, int subsample_, int ntargets, float * xpos,
                               float * ypos, float * zlambda, float * mag, float zerop, long * sizes,
                               int Npts, int device)

            void publish()


IF USE_MPI == 1:
    # from mpi4py import MPI
    from mpi4py cimport MPI
    # C-level cdef, typed, Python objects
    # from mpi4py cimport mpi_c as mpi
    from mpi4py cimport libmpi as mpi

IF USE_MPI == 2:
    cimport mpi4py.MPI as MPI
    cimport mpi4py.libmpi as mpi


IF USE_MPI:
    cdef inline Bcast(carma_obj[float] * obj, int root):
        """Broadcast the content of a carma_obj<float>

        :parameters:
            obj: (carma_obj<float>) : carma_obj to broadcast

            root: (int) : root of the MPI broadcast
        """
        cdef int i
        cdef int size = < int > obj.getNbElem()

        cdef float * ptr
        ptr = < float * > malloc(size * sizeof(float))

        obj.device2host(ptr)

        mpi.MPI_Bcast(ptr, size, mpi.MPI_FLOAT, root, mpi.MPI_COMM_WORLD)

        obj.host2device(ptr)

        free(ptr)

    cdef inline Bcast_cudaAware(carma_obj[float] * obj, int root):
        """Broadcast the content of a carma_obj<float>
           Using cuda_aware

        :parameters:
            obj: (carma_obj<float>) : carma_obj to broadcast

            root: (int) : root of the MPI broadcast
        """

        cdef int i
        cdef int size = < int > obj.getNbElem()

        cdef float * ptr
        ptr = obj.getData()

        mpi.MPI_Bcast(ptr, size, mpi.MPI_FLOAT, root, mpi.MPI_COMM_WORLD)
