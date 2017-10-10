from cython.operator cimport dereference as deref, preincrement as inc

import numpy as np

import pandas as pd
from scipy import interpolate
import shesha_constants as scons

# import astropy.io.fits as pfits
# max_extent signature
#################################################
# P-Class Dms
#################################################
cdef class Dms:
    def __cinit__(self, naga_context context, long ndm):
        self.context = context
        self.dms = new sutra_dms(ndm)

    def __dealloc__(self):
        del self.dms

    def __str__(self):
        info = "DMs object:\n"
        info += "Contains " + str(self.dms.d_dms.size()) + " DMs:\n"
        info += "DM # | Type  |   Alt   | Nact | Dim\n"
        cdef vector[sutra_dm * ].iterator it_dms = self.dms.d_dms.begin()
        cdef vector[type_screen].iterator it_type = self.dms.d_type.begin()
        cdef sutra_dm * dm
        cdef type_screen ts
        cdef int i = 1
        while(it_dms != self.dms.d_dms.end()):
            dm = deref(it_dms)
            ts = deref(it_type)
            info += "%4d" % i + " | " + "%5s" % ts.first + " | " + "%7d" % ts.second + \
                " | " + "%4d" % dm.ninflu + " | " + "%4d" % dm.dim + "\n"
            i = i + 1
            inc(it_dms)
            inc(it_type)
        info += "--------------------------------------------------------"
        return info

    def add_dm(self, bytes type, float alt, long dim, long ninflu, long influsize, long ninflupos, long npts, float push4imat, int device=-1):
        """Add a dm into a Dms object

        :parameters:
            type: (str) : dm type to remove,

            alt: (float) : dm conjugaison altitude to remove,

            ninflu: (long) : ,

            influsize: (long) : ,

            ninflupos: (long) : ,

            npts: (long) : ,

            push4imat: (float) : ,

            device: (int) : device where the DM will be create (default=-1):

        """
        if(device > -1):
            self.context.set_activeDevice(device, 1)
        else:
            device = self.context.get_activeDevice()

        self.dms.add_dm(self.context.c, type, alt, dim, ninflu,
                        influsize, ninflupos, npts, push4imat, device)

    def remove_dm(self, bytes type, float alt):
        """Remove a dm from a Dms object

        :parameters:
            type: (str) : dm type to remove

            alt: (float) : dm conjugaison altitude to remove
        """
        self.dms.remove_dm(type, alt)

    def resetdm(self, bytes type, float alt):
        """Reset the shape of the DM to 0

        :parameters:
            type: (str) : dm type

            alt: (float) : dm conjugaison altitude
        """
        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            raise StandardError("Error in reset dm ")

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        self.dms.d_dms[inddm].reset_shape()

    def oneactu(self, bytes type, float alt, int nactu, float ampli):
        """Push on on the nactu actuator of the DM with ampli amplitude and compute
        the corresponding shape

        :parameters:
            type: (str) : dm type

            alt: (float) : dm conjugaison altitude

            nactu: (int) : actuator number

            ampli: (float): amplitude
        """
        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            raise StandardError("One actuator error")

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        self.dms.d_dms[inddm].comp_oneactu(nactu, ampli)

    def load_pzt(self, float alt,
                 np.ndarray[ndim=3, dtype=np.float32_t] influ,
                 np.ndarray[ndim=1, dtype=np.int32_t] influpos,
                 np.ndarray[ndim=1, dtype=np.int32_t] npoints,
                 np.ndarray[ndim=1, dtype=np.int32_t] istart,
                 np.ndarray[ndim=1, dtype=np.int32_t] xoff,
                 np.ndarray[ndim=1, dtype=np.int32_t] yoff):
        """Load all the arrays computed during the initialization
        for a pzt DM in a sutra_dms object

        :parameters:
            alt: (float) : dm conjugaison altitude

            influ: (np.ndarray[ndim=3,dtype=np.float32_t]) : influence functions

            influpos: (np.ndarray[ndim=1,dtype=np.int32_t]) : positions of the IF

            npoints: (np.ndarray[ndim=1,dtype=np.int32_t]) : for each pixel on the DM screen,
                                                            the number of IF which impact on this pixel

            istart: (np.ndarray[ndim=1,dtype=np.int32_t]) :

            xoff: (np.ndarray[ndim=1,dtype=np.int32_t]) : x-offset

            yoff: (np.ndarray[ndim=1,dtype=np.int32_t]) :y-offset

        """

        cdef np.ndarray[ndim = 3, dtype = np.float32_t] influ_F = influ.T.copy()
        cdef np.ndarray[dtype = np.int32_t] npoints_F = npoints.T.copy()

        cdef int inddm = self.dms.get_inddm(scons.DmType.PZT, alt)
        if(inddm < 0):
            err = "unknown error whith load_pzt\nDM (pzt" + str(
                alt) + ") doesn't exist"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)

        cdef int ntot_influpos = influpos.shape[0]
        cdef int ntot_istart = istart.shape[0]
        cdef int n = self.dms.d_dms[inddm].influsize * \
            self.dms.d_dms[inddm].influsize

        cdef float * influ2
        influ2 = < float * > malloc(ntot_influpos * sizeof(float))
        cdef tuple_t[float] * influ3 = NULL
        cdef int * influpos2
        influpos2 = < int * > malloc(ntot_influpos * sizeof(int))
        cdef int * istart2
        istart2 = < int * > malloc((ntot_istart + 1) * sizeof(int))

        cdef int i
        for i in range(ntot_influpos):
            influ2[i] = influpos[i] / n

        for i in range(ntot_istart):
            istart2[i] = istart[i]
        istart2[ntot_istart] = istart2[
            ntot_istart - 1] + npoints[ntot_istart - 1]

        # TODO preprocessing: COMPN==2 or 3
        """
#if(COMPN == 2)
      influ3 = new struct tuple_t<float>[ntot_influpos];

      for(int  i = 0; i < ntot_influpos; i++)
    	  influ3[i] = {influpos2[i], influ2[i]};

#elif(COMPN == 3)
      influ3 = new struct tuple_t<float>[ntot_istart * MAXSPOT];

      //For each pixel of screen
      int count = 0;
      for(int  i = 0; i < ntot_istart; i++)
      {
    	  int j = 0;
    	  //For each influence function, cpy the value of postition and influ
    	  for(; j < npoints[i]; j++){
    		  influ3[i * MAXSPOT + j] = {influpos2[count], influ2[count]};
    		  count++;
    	  }
    	  //Fill the last element with 0
    	  for(; j < MAXSPOT; j++)
    		  influ3[i * MAXSPOT + j] = {0, 0.0};
      }
#endif
        """

        self.dms.d_dms[inddm].pzt_loadarrays( < float * > influ_F.data,
                                             influ2,
                                             influ3,
                                             < int * > influpos.data,
                                             influpos2,
                                             < int * > npoints_F.data,
                                             istart2,
                                             < int * > xoff.data,
                                             < int * > yoff.data)
        free(influ2)
        free(influpos2)
        free(istart2)

    def load_kl(self, float alt,
                np.ndarray[ndim=2, dtype=np.float32_t] rabas,
                np.ndarray[ndim=2, dtype=np.float32_t] azbas,
                np.ndarray[ndim=1, dtype=np.int32_t] ord,
                np.ndarray[ndim=2, dtype=np.float32_t] cr,
                np.ndarray[ndim=2, dtype=np.float32_t] cp):
        """Load all the arrays computed during the initialization
        for a kl DM in a sutra_dms object

        :parameters:
            alt: (float) : dm conjugaison altitude

            rabas: (np.ndarray[ndim=1,dtype=np.float32_t]) : TODO

            azbas: (np.ndarray[ndim=1,dtype=np.float32_t]) :

            ord: (np.ndarray[ndim=1,dtype=np.int32_t]) :

            cr: (np.ndarray[ndim=1,dtype=np.float32_t]) :

            cp: (np.ndarray[ndim=1,dtype=np.float32_t]) :
        """

        cdef int inddm = self.dms.get_inddm(scons.DmType.KL, alt)
        if(inddm < 0):
            err = "unknown error whith load_kl\nDM (kl" + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)

        cdef np.ndarray[ndim = 2, dtype = np.float32_t] rabas_F = rabas.T.copy()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] cr_F = cr.T.copy()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] cp_F = cp.T.copy()

        self.dms.d_dms[inddm].kl_loadarrays( < float * > rabas_F.data,
                                            < float * > azbas.data,
                                            < int * > ord.data,
                                            < float * > cr_F.data,
                                            < float * > cp_F.data)

    def load_tt(self, float alt, np.ndarray[ndim=3, dtype=np.float32_t] influ):
        """Load all the arrays computed during the initialization
        for a tt DM in a sutra_dms object

        :parameters:
            alt: (float) : dm conjugaison altitude

            influ: (np.ndarray[ndim=3,dtype=np.float32_t]) : influence functions
        """
        cdef int inddm = self.dms.get_inddm(scons.DmType.TT, alt)
        if(inddm < 0):
            err = "unknown error whith load_tt\nDM (tt" + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] influ_F = influ.T.copy()
        self.dms.d_dms[inddm].d_influ.host2device( < float * > influ_F.data)

    def shape_dm(self, bytes type, float alt):
        """Compute the shape of the DM in a sutra_dm object

            type: (str) : dm type

            alt: (float) : dm conjugaison altitude
        """
        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            err = "unknown error whith set_comm\nDM (" + type + ", " + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)

        self.dms.d_dms[inddm].comp_shape()

    def compute_KLbasis(self, bytes type, float alt,
                        np.ndarray[ndim=1, dtype=np.float32_t] xpos, np.ndarray[ndim=1, dtype=np.float32_t] ypos,
                        np.ndarray[ndim=1, dtype=np.int32_t] indx_pup, long dim, float norm, float ampli):
        """Compute a Karhunen-Loeve basis for the dm:
            - compute the phase covariance matrix on the actuators using Kolmogorov
            - compute the geometric covariance matrix
            - double diagonalisation to obtain KL basis

        :parameters:
            type: (str) : dm type

            alt: (float) : dm conjugaison altitude

            xpos: (np.ndarray[ndim=1,dtype=np.float32_t]) : x-position of actuators

            ypos: (np.ndarray[ndim=1,dtype=np.float32_t]) : y-position of actuators

            indx_pup: (np.ndarray[ndim=1,dtype=np.int32_t]) : indices of where(pup)

            dim: (long) : number of where(pup)

            norm: (float) : normalization factor

            ampli: (float) : amplitude
        """

        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            err = "unknown error with computeKLbasis function\nDM(" + type + "," + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)

        self.dms.d_dms[inddm].compute_KLbasis(< float * > xpos.data, < float * > ypos.data,
                                               < int * > indx_pup.data, dim, norm, ampli)

    def comp_oneactu(self, bytes type, float alt, int nactu, float ampli):
        """Compute the shape of the dm when pushing the nactu actuator

        :parameters:
            type: (str) : dm type

            alt: (float) : dm conjugaison altitude

            nactu: (int) : actuator number pushed

            ampli: (float) : amplitude
        """

        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            err = "unknown error with get_dm function DM(" + \
                type + "," + alt + ") doesnt exists"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        self.dms.d_dms[inddm].comp_oneactu(nactu, ampli)

    """
           _______. _______ .___________.                _______  _______ .___________.
          /       ||   ____||           |     ___       /  _____||   ____||           |
         |   (----`|  |__   `---|  |----`    ( _ )     |  |  __  |  |__   `---|  |----`
          \   \    |   __|      |  |         / _ \/\   |  | |_ | |   __|      |  |
      .----)   |   |  |____     |  |        | (_>  <   |  |__| | |  |____     |  |
      |_______/    |_______|    |__|         \___/\/    \______| |_______|    |__|

    """

    def set_full_comm(self, np.ndarray[ndim=1, dtype=np.float32_t] comm,
                      bool shape_dm=True):
        """Set the voltage command

            comm: (np.ndarray[ndim=1,dtype=np.float32_t]) : voltage vector

            shape_dm: (bool) : perform the dm_shape after the load (default=True)
        """
        nact_total = self.dms.nact_total()
        if nact_total != comm.size:
            raise ValueError("Incorrect size of voltage vector")

        cdef int comm_index = 0
        cdef float * comm_data = < float * > comm.data
        for inddm in range(self.dms.ndm()):
            self.dms.d_dms[inddm].d_comm.host2device(& comm_data[comm_index])
            comm_index += self.dms.d_dms[inddm].nact()
            if shape_dm:
                self.dms.d_dms[inddm].comp_shape()

    def set_comm(self, bytes type, float alt,
                 np.ndarray[ndim=1, dtype=np.float32_t] comm,
                 bool shape_dm=False):
        """Set the voltage command on a sutra_dm

            type: (str) : dm type

            alt: (float) : dm conjugaison altitude

            comm: (np.ndarray[ndim=1,dtype=np.float32_t]) : voltage vector

            shape_dm: (bool) : perform the dm_shape after the load (default=False)
        """

        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            err = "unknown error whith set_comm\nDM (" + type + ", " + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        self.dms.d_dms[inddm].d_comm.host2device( < float * > comm.data)
        if shape_dm:
            self.dms.d_dms[inddm].comp_shape()
        if shape_dm:
            self.dms.d_dms[inddm].comp_shape()

    def get_KLbasis(self, bytes type, float alt):
        """Return the klbasis computed by computeKLbasis

        :parameters:
            type: (str) : dm type

            alt: (float) : dm conjugaison altitude
        :return:
            KLbasis : (np.ndarray(dims=2,dtype=np.float32)) : the KL basis
        """

        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            err = "unknown error with get_KLbasis function DM(" + \
                type + "," + alt + ") doesnt exists"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        cdef const long * dims = self.dms.d_dms[inddm].d_KLbasis.getDims()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)

        self.dms.d_dms[inddm].d_KLbasis.device2host( < float * > data_F.data)
        return data_F.T.copy()

    def get_dm(self, bytes type, float alt):
        """Return the shape of the dm

        :parameters:
            type: (str) : dm type

            alt: (float) : dm conjugaison altitude
        :return:
            data : (np.ndarray(dims=2,dtype=np.float32)) : DM shape
        """

        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            err = "unknown error with get_dm function DM(" + \
                type + "," + alt + ") doesnt exists"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        cdef const long * dims = self.dms.d_dms[inddm].d_shape.d_screen.getDims()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)

        self.dms.d_dms[inddm].d_shape.d_screen.device2host( < float * > data_F.data)
        return data_F.T.copy()

    def get_comm(self, bytes type, float alt):
        """Return the voltage command of the sutra_dm

        :parameters:
            type: (str) : dm type

            alt: (float) : dm conjugaison altitude
g        :return:
            data : (np.ndarray(dims=1,dtype=np.float32)) : voltage vector
        """

        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            err = "unknown error with get_dm function DM(" + \
                type + "," + alt + ") doesnt exists"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        cdef const long * dims = self.dms.d_dms[inddm].d_comm.getDims()
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data = np.zeros((dims[1]), dtype=np.float32)

        self.dms.d_dms[inddm].d_comm.device2host( < float * > data.data)
        return data

    def get_influ(self, bytes type, float alt):
        """Return the influence functions of the DM

        :parameters:
            type: (str) : dm type

            alt: (float) : dm conjugaison altitude
        :return:
            data : (np.ndarray(dims=3,dtype=np.float32)) : influence functions
        """

        cdef int inddm = self.dms.get_inddm(type, alt)
        if(inddm < 0):
            err = "unknown error with get_dm function DM(" + \
                type + "," + alt + ") doesnt exists"
            raise ValueError(err)

        self.context.set_activeDevice(self.dms.d_dms[inddm].device, 1)

        cdef const long * dims = self.dms.d_dms[inddm].d_influ.getDims()
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] data_F = np.zeros((dims[3], dims[2], dims[1]), dtype=np.float32)

        self.dms.d_dms[inddm].d_influ.device2host( < float * > data_F.data)
        return data_F.T.copy()

    def get_IFsparse(self, bytes type, float alt, np.ndarray[ndim=1, dtype=np.int32_t] indx_pup):
        """Returns the influence functions matrix of a pzt DM as a sparse matrix
        Return a scipy.sparse object which shape is (nactus,Npts in the pupil)

        :parameters:
            type: (str) : DM type
            alt: (float) : DM altitude
            indx_pup: (np.ndarray[ndim=1, dtype=np.int32_t]) : valid indices of the pupil in the DM support
        :return:
            IF : (scipy.sparse) : influence functions matrix
        """

        cdef int inddm = self.dms.get_inddm(type, alt)
        self.context.set_activeDeviceForCpy(self.dms.d_dms[inddm].device, 1)

        cdef carma_sparse_obj[double] * d_IFsparse = NULL
        cdef carma_obj[int] * d_indx
        cdef long dims[2]
        dims[0] = 1
        dims[1] = indx_pup.size

        if(type == scons.DmType.PZT):
            d_indx = new carma_obj[int](self.context.c, dims, < int * > indx_pup.data)
            sparse = naga_sparse_obj_Double()
            self.dms.d_dms[inddm].get_IF_sparse(
                d_IFsparse, d_indx.getData(), indx_pup.size, float(1.0), 1)
            sparse.copy(d_IFsparse)
            del d_indx
            del d_IFsparse
            return sparse.get_sparse()
        else:
            raise ValueError(
                "This function only works with pzt DM (tt influence functions are not sparse)")

    def get_IFtt(self, bytes type, float alt, np.ndarray[ndim=1, dtype=np.int32_t] indx_pup):
        """Returns the influence functions matrix of a tt DM

        :parameters:
            type: (str) : DM type
            alt: (float) : DM altitude
            indx_pup: (np.ndarray[ndim=1, dtype=np.int32_t]) : valid indices of the pupil in the DM support
        :return:
            IFtt : (np.ndarray[ndim=2, dtype=np.float32_t]) : influence functions matrix
        """

        cdef int inddm = self.dms.get_inddm(type, alt)
        self.context.set_activeDeviceForCpy(self.dms.d_dms[inddm].device, 1)

        cdef carma_obj[int] * d_indx
        cdef carma_obj[float] * d_IFtt
        cdef long dims[2]
        dims[0] = 1
        dims[1] = indx_pup.size
        cdef long dims2[3]
        dims2[0] = 2
        dims2[1] = indx_pup.size
        dims2[2] = 2
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F

        if(type == scons.DmType.TT):
            d_indx = new carma_obj[int](self.context.c, dims, < int * > indx_pup.data)
            d_IFtt = new carma_obj[float](self.context.c, dims2)
            self.dms.d_dms[inddm].get_IF(
                d_IFtt.getData(),
                d_indx.getData(),
                indx_pup.size,
                float(1.0))
            data_F = np.zeros((dims2[2], dims2[1]), dtype=np.float32)
            d_IFtt.device2host( < float * > data_F.data)

            del d_indx
            del d_IFtt
            return data_F.T.copy()
        else:
            raise ValueError(
                "This function only works with tt DM (for pzt, use get_IFsparse)")
