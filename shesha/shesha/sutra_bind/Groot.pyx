import numpy as np

from naga.context import context as naga_context


def groot_init(int nactus, int nlayers, float gsangle, float fc,
               np.ndarray[ndim=1, dtype=np.float32_t] vdt,
               np.ndarray[ndim=1, dtype=np.float32_t] Htheta,
               np.ndarray[ndim=1, dtype=np.float32_t] L0,
               np.ndarray[ndim=1, dtype=np.float32_t] winddir,
               np.ndarray[ndim=1, dtype=np.float32_t] scale,
               np.ndarray[ndim=1, dtype=np.float32_t] xpos,
               np.ndarray[ndim=1, dtype=np.float32_t] ypos,
               np.ndarray[ndim=2, dtype=np.float32_t] pzt2tt,
               np.ndarray[ndim=2, dtype=np.float32_t] TTPfilter,
               np.ndarray[ndim=2, dtype=np.float32_t] Nact):
    """ Initialize a Groot object for residual covariance matrix estimation
    :parameters:
        nactus : (int) : number of actuators
        nlayers : (int) : number of turbulent layers
        gsangle : (float) : Guide star angle (np.arctan2(wfs.ypos,wfs.xpos)) [rad]
        vdt : (np.ndarray[ndim=1,dtype=float32_t]) : windspeed * iter_time / gain for each layer
        Htheta : (np.ndarray[ndim=1,dtype=float32_t]) : Distance of the guide star from the axis at altitude h [m]
        L0 : (np.ndarray[ndim=1,dtype=float32_t]) : Outer scale of each layers [m]
        winddir : (np.ndarray[ndim=1,dtype=float32_t]) : Wind directions for each layer
        scale : (np.ndarray[ndim=1,dtype=float32_t]) : Scaling factor (r0**(-5/3) * frac * (lambda/2pi)**2
        xpos : (np.ndarray[ndim=1,dtype=float32_t]) : X-positions of pzt DM actuators [m]
        ypos : (np.ndarray[ndim=1,dtype=float32_t]) : Y-positions of pzt DM actuators [m]
        pzt2tt : (np.ndarray[ndim=2,dtype=float32_t]) : pzt to tt actuators projection matrix
        TTPfilter : (np.ndarray[ndim=2,dtype=float32_t]) : tip-tilt+piston filter matrix (Btt.dot(P))
        Nact : (np.ndarray[ndim=2,dtype=float32_t]) : pzt actuators coupling matrix
    """

    return Groot(nactus, nlayers, gsangle, fc, vdt, Htheta, L0, winddir, scale,
                 xpos, ypos, pzt2tt, TTPfilter, Nact)


def groot_init_alias(int nvalid, float fc, float scale, float d, int npts,
                     np.ndarray[ndim=1, dtype=np.float32_t] xpos,
                     np.ndarray[ndim=1, dtype=np.float32_t] ypos,
                     np.ndarray[ndim=1, dtype=np.float32_t] weights):
    """ Initializes a Groot object dedicated to aliasing modeling
    """

    return GrootAlias(nvalid, fc, scale, d, npts, xpos, ypos, weights)


cdef class GrootAlias:

    def __cinit__(self, int nvalid, float fc, float scale, float d, int npts,
                  np.ndarray[ndim=1, dtype=np.float32_t] xpos,
                  np.ndarray[ndim=1, dtype=np.float32_t] ypos,
                  np.ndarray[ndim=1, dtype=np.float32_t] weights):
        """ Initialize a shesha_groot object for aliasing modeling
        """
        cdef carma_context * context = &carma_context.instance()
        cdef int device
        device = context.get_activeDevice()
        self.device = device

        self.groot = new sutra_groot(context, device, nvalid, < float * > weights.data, scale, < float * > xpos.data, < float * > ypos.data, fc, d, npts)

    def __dealloc__(self):
        del self.groot

    def compute_Calias(self):
        """ Compute the aliasing covariance matrix with GROOT model
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.groot.device, 1)

        self.groot.compute_Calias()

    def get_Calias(self):
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.groot.device, 1)
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.groot.d_CaXX.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.groot.d_CaXX.device2host(< float * > data_F.data)
        data = np.zeros((2 * dims[1], 2 * dims[1]), dtype=np.float32)
        data[:dims[1], :dims[1]] = np.reshape(
            data_F.flatten("F"), (dims[1], dims[2]))

        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.groot.d_CaYY.device2host(< float * > data_F.data)
        data[dims[1]:, dims[1]:] = np.reshape(
            data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_CaXX(self):
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.groot.device, 1)
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.groot.d_CaXX.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.groot.d_CaXX.device2host(< float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_CaYY(self):
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.groot.device, 1)
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.groot.d_CaYY.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.groot.d_CaYY.device2host(< float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

cdef class Groot:

    def __cinit__(self, int nactus, int nlayers, float gsangle, float fc,
                  np.ndarray[ndim=1, dtype=np.float32_t] vdt,
                  np.ndarray[ndim=1, dtype=np.float32_t] Htheta,
                  np.ndarray[ndim=1, dtype=np.float32_t] L0,
                  np.ndarray[ndim=1, dtype=np.float32_t] winddir,
                  np.ndarray[ndim=1, dtype=np.float32_t] scale,
                  np.ndarray[ndim=1, dtype=np.float32_t] xpos,
                  np.ndarray[ndim=1, dtype=np.float32_t] ypos,
                  np.ndarray[ndim=2, dtype=np.float32_t] pzt2tt,
                  np.ndarray[ndim=2, dtype=np.float32_t] TTPfilter,
                  np.ndarray[ndim=2, dtype=np.float32_t] Nact):
        """ Initialize a shesha_groot object for residual covariance matrix estimation
        :parameters:
          nactus : (int) : number of actuators
          nlayers : (int) : number of turbulent layers
          gsangle : (float) : Guide star angle (np.arctan2(wfs.ypos,wfs.xpos)) [rad]
          vdt : (np.ndarray[ndim=1,dtype=float32_t]) : windspeed * iter_time / gain for each layer
          Htheta : (np.ndarray[ndim=1,dtype=float32_t]) : Distance of the guide star from the axis at altitude h [m]
          L0 : (np.ndarray[ndim=1,dtype=float32_t]) : Outer scale of each layers [m]
          winddir : (np.ndarray[ndim=1,dtype=float32_t]) : Wind directions for each layer
          scale : (np.ndarray[ndim=1,dtype=float32_t]) : Scaling factor (r0**(-5/3) * frac * (lambda/2pi)**2
          xpos : (np.ndarray[ndim=1,dtype=float32_t]) : X-positions of pzt DM actuators [m]
          ypos : (np.ndarray[ndim=1,dtype=float32_t]) : Y-positions of pzt DM actuators [m]
          pzt2tt : (np.ndarray[ndim=2,dtype=float32_t]) : pzt to tt actuators projection matrix
          TTPfilter : (np.ndarray[ndim=2,dtype=float32_t]) : tip-tilt+piston filter matrix (Btt.dot(P))
          Nact : (np.ndarray[ndim=2,dtype=float32_t]) : pzt actuators coupling matrix
        """

        cdef carma_context * context = &carma_context.instance()
        cdef int device
        device = context.get_activeDevice()
        self.device = device
        cdef np.ndarray[dtype= np.float32_t] pzt2tt_F = pzt2tt.flatten("F")
        cdef np.ndarray[dtype= np.float32_t] TTPfilter_F = TTPfilter.flatten("F")
        cdef np.ndarray[dtype= np.float32_t] Nact_F = Nact.flatten("F")

        self.groot = new sutra_groot(context, device, nactus, nlayers,
                                     gsangle, < float * > vdt.data,
                                     < float * > Htheta.data, < float * > L0.data,
                                     < float * > winddir.data, < float * > scale.data,
                                     < float * > pzt2tt_F.data, < float * > TTPfilter_F.data,
                                     < float * > Nact_F.data, < float * > xpos.data,
                                     < float * > ypos.data, fc)

    def __dealloc__(self):
        del self.groot

    def compute_Cerr(self):
        """ Compute the residual error covariance matrix Cerr with GROOT model
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.groot.device, 1)

        self.groot.compute_Cerr()

    def get_Cerr(self):
        """Return Cerr computed on the pzt DM (nactus x nactus).
        Use get_TTcomp to retrieve the TT componenent and merge the two matrix
        to have the full Cerr

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : Cerr
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.groot.device, 1)
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.groot.d_Cerr.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.groot.d_Cerr.device2host(< float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_TTcomp(self):
        """Return Cerr computed on the TT DM (2 x 2).
        Use get_Cerr to retrieve the pzt componenent and merge the two matrix
        to have the full Cerr

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : TTcomp
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.groot.device, 1)
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.groot.d_TT.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.groot.d_TT.device2host(< float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_tab_int_x(self):
        """Return the X tabulation for integral

        :return:
            data : (np.ndarray[ndim=1,dtype=np.int32_t]) : tab_int_x
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.groot.device, 1)
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.groot.d_tab_int_x.getDims()
        data = np.zeros(dims[1], dtype=np.float32)
        self.groot.d_tab_int_x.device2host(< float * > data.data)

        return data

    def get_tab_int_y(self):
        """Return the Y tabulation for integral

        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : tab_int_y
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.groot.device, 1)
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.groot.d_tab_int_y.getDims()
        data = np.zeros(dims[1], dtype=np.float32)
        self.groot.d_tab_int_y.device2host(< float * > data.data)

        return data
