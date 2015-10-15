

"""
cdef extern from "carma_obj.h":
    int carma_svd[T](carma_obj[T] *imat, carma_obj[T] *eigenvals,
        carma_obj[T] *mod2act, carma_obj[T] *mes2mod)
    int carma_getri[T](carma_obj[T] *d_iA)
    int carma_potri[T](carma_obj[T] *d_iA)

    int carma_syevd[T](char jobz, long N, T *mat, T *eigenvals)

cdef extern from "carma_host_obj.h":
    int carma_svd_cpu[T](carma_host_obj[T] *imat,
        carma_host_obj[T] *eigenvals, carma_host_obj[T] *mod2act,
        carma_host_obj[T] *mes2mod)
    int carma_getri_cpu[T](long N, T *h_A)
    int carma_potri_cpu[T](long N, T *h_A)
    int carma_syevd_cpu[T](char jobz, int N, T *h_A, T *eigenvals)
"""


def svd_host_Double(chakra_host_obj_Double2D mat,
                    chakra_host_obj_Double1D eigenvals,
                    chakra_host_obj_Double2D U,
                    chakra_host_obj_Double2D VT):
    """Call carma_svd_cpu
    
    chakra_host_obj_Double2D mat:
    chakra_host_obj_Double1D eigenvals:
    chakra_host_obj_Double2D U:
    chakra_host_obj_Double2D VT:
    """
    carma_svd_cpu[double](mat.c_h, eigenvals.c_h, VT.c_h, U.c_h)

def svd_Double(chakra_obj_Double2D mat,
                chakra_obj_Double1D eigenvals,
                chakra_obj_Double2D U,
                chakra_obj_Double2D VT):
    """Call carma_svd
    """

    carma_svd[double](mat.c_o, eigenvals.c_o, VT.c_o, U.c_o)


def getri_Double(chakra_obj_Double2D d_mat):
    return carma_getri[double](d_mat.c_o)


def getri_host_Double(chakra_host_obj_Double2D h_mat):
    cdef const long *dims=h_mat.c_h.getDims()
    return carma_getri_cpu[double](dims[1],h_mat.data_h)


def potri_Double(chakra_obj_Double2D d_mat ):
    return carma_potri[double](d_mat.c_o)


def potri_host_Double(chakra_host_obj_Double2D h_mat):
    cdef const long *dims=h_mat.c_h.getDims()
    return carma_potri_cpu[double](dims[1],h_mat.c_h.getData())


def syevd_Double(chakra_obj_Double2D d_A,
                np.ndarray[ndim=1,dtype=np.float64_t] eigenvals,
                chakra_obj_Double2D U=None,
                bool computeU=True):

    cdef long N=eigenvals.shape[0]
    
    if(U is None):
        if(computeU==True):
            carma_syevd[double]('V',N,d_A.c_o.getData(), <double*>eigenvals.data)
        else:
            carma_syevd[double]('N',N,d_A.c_o.getData(), <double*>eigenvals.data)
    else:
        U.c_o.copy(d_A.c_o,1,1)
        if(computeU==True):
            carma_syevd[double]('V',N,U.c_o.getData(), <double*>eigenvals.data)
        else:
            carma_syevd[double]('N',N,U.c_o.getData(), <double*>eigenvals.data)



def syevd_host_Double(chakra_host_obj_Double2D h_A,
                np.ndarray[ndim=1,dtype=np.float64_t] eigenvals,
                chakra_host_obj_Double2D U=None,
                bool computeU=True):
    
    cdef long N=eigenvals.shape[0]

    if(U is None):
        if(computeU==True):
            carma_syevd_cpu[double]('V',N, h_A.data_h, <double*>eigenvals.data)
        else:
            carma_syevd_cpu[double]('N',N, h_A.data_h, <double*>eigenvals.data)
    else:
        U.fill_from(h_A)
        if(computeU==True):
            carma_syevd_cpu[double]('V',N,U.data_h, <double*>eigenvals.data)
        else:
            carma_syevd_cpu[double]('N',N,U.data_h, <double*>eigenvals.data)




def svd_host_Double(chakra_host_obj_Double2D mat,
                    chakra_host_obj_Double1D eigenvals,
                    chakra_host_obj_Double2D U,
                    chakra_host_obj_Double2D VT):
    """Call carma_svd_cpu
    
    chakra_host_obj_Double2D mat:
    chakra_host_obj_Double1D eigenvals:
    chakra_host_obj_Double2D U:
    chakra_host_obj_Double2D VT:
    """
    carma_svd_cpu[double](mat.c_h, eigenvals.c_h, VT.c_h, U.c_h)

def svd_Double(chakra_obj_Double2D mat,
                chakra_obj_Double1D eigenvals,
                chakra_obj_Double2D U,
                chakra_obj_Double2D VT):
    """Call carma_svd
    """

    carma_svd[double](mat.c_o, eigenvals.c_o, VT.c_o, U.c_o)


def getri_Double(chakra_obj_Double2D d_mat):
    return carma_getri[double](d_mat.c_o)


def getri_host_Double(chakra_host_obj_Double2D h_mat):
    cdef const long *dims=h_mat.c_h.getDims()
    return carma_getri_cpu[double](dims[1],h_mat.data_h)


def potri_Double(chakra_obj_Double2D d_mat ):
    return carma_potri[double](d_mat.c_o)


def potri_host_Double(chakra_host_obj_Double1D h_mat):
    cdef const long *dims=h_mat.c_h.getDims()
    return carma_potri_cpu[double](dims[1],h_mat.c_h.getData())


def syevd_Double(chakra_obj_Double2D d_A,
                np.ndarray[ndim=1,dtype=np.float64_t] eigenvals,
                chakra_obj_Double2D U=None,
                bool computeU=True):

    cdef long N=eigenvals.shape[0]
    
    if(U is None):
        if(computeU==True):
            carma_syevd[double]('V',N,d_A.c_o.getData(), <double*>eigenvals.data)
        else:
            carma_syevd[double]('N',N,d_A.c_o.getData(), <double*>eigenvals.data)
    else:
        U.c_o.copy(d_A.c_o,1,1)
        if(computeU==True):
            carma_syevd[double]('V',N,U.c_o.getData(), <double*>eigenvals.data)
        else:
            carma_syevd[double]('N',N,U.c_o.getData(), <double*>eigenvals.data)



def syevd_host_Double(chakra_host_obj_Double2D h_A,
                np.ndarray[ndim=1,dtype=np.float64_t] eigenvals,
                chakra_host_obj_Double2D U=None,
                bool computeU=True):
    
    cdef long N=eigenvals.shape[0]

    if(U is None):
        if(computeU==True):
            carma_syevd_cpu[double]('V',N, h_A.data_h, <double*>eigenvals.data)
        else:
            carma_syevd_cpu[double]('N',N, h_A.data_h, <double*>eigenvals.data)
    else:
        U.fill_from(h_A)
        if(computeU==True):
            carma_syevd_cpu[double]('V',N,U.data_h, <double*>eigenvals.data)
        else:
            carma_syevd_cpu[double]('N',N,U.data_h, <double*>eigenvals.data)




def svd_host_Float(chakra_host_obj_Float2D mat,
                    chakra_host_obj_Float1D eigenvals,
                    chakra_host_obj_Float2D U,
                    chakra_host_obj_Float2D VT):
    """Call carma_svd_cpu
    
    chakra_host_obj_Float2D mat:
    chakra_host_obj_Float1D eigenvals:
    chakra_host_obj_Float2D U:
    chakra_host_obj_Float2D VT:
    """
    carma_svd_cpu[float](mat.c_h, eigenvals.c_h, VT.c_h, U.c_h)

def svd_Float(chakra_obj_Float2D mat,
                chakra_obj_Float1D eigenvals,
                chakra_obj_Float2D U,
                chakra_obj_Float2D VT):
    """Call carma_svd
    """

    carma_svd[float](mat.c_o, eigenvals.c_o, VT.c_o, U.c_o)


def getri_Float(chakra_obj_Float2D d_mat):
    return carma_getri[float](d_mat.c_o)


def getri_host_Float(chakra_host_obj_Float2D h_mat):
    cdef const long *dims=h_mat.c_h.getDims()
    return carma_getri_cpu[float](dims[1],h_mat.data_h)


def potri_Float(chakra_obj_Float2D d_mat ):
    return carma_potri[float](d_mat.c_o)


def potri_host_Float(chakra_host_obj_Float1D h_mat):
    cdef const long *dims=h_mat.c_h.getDims()
    return carma_potri_cpu[float](dims[1],h_mat.c_h.getData())


def syevd_Float(chakra_obj_Float2D d_A,
                np.ndarray[ndim=1,dtype=np.float32_t] eigenvals,
                chakra_obj_Float2D U=None,
                bool computeU=True):

    cdef long N=eigenvals.shape[0]
    
    if(U is None):
        if(computeU==True):
            carma_syevd[float]('V',N,d_A.c_o.getData(), <float*>eigenvals.data)
        else:
            carma_syevd[float]('N',N,d_A.c_o.getData(), <float*>eigenvals.data)
    else:
        U.c_o.copy(d_A.c_o,1,1)
        if(computeU==True):
            carma_syevd[float]('V',N,U.c_o.getData(), <float*>eigenvals.data)
        else:
            carma_syevd[float]('N',N,U.c_o.getData(), <float*>eigenvals.data)



def syevd_host_Float(chakra_host_obj_Float2D h_A,
                np.ndarray[ndim=1,dtype=np.float32_t] eigenvals,
                chakra_host_obj_Float2D U=None,
                bool computeU=True):
    
    cdef long N=eigenvals.shape[0]

    if(U is None):
        if(computeU==True):
            carma_syevd_cpu[float]('V',N, h_A.data_h, <float*>eigenvals.data)
        else:
            carma_syevd_cpu[float]('N',N, h_A.data_h, <float*>eigenvals.data)
    else:
        U.fill_from(h_A)
        if(computeU==True):
            carma_syevd_cpu[float]('V',N,U.data_h, <float*>eigenvals.data)
        else:
            carma_syevd_cpu[float]('N',N,U.data_h, <float*>eigenvals.data)




def svd_host_Float(chakra_host_obj_Float2D mat,
                    chakra_host_obj_Float1D eigenvals,
                    chakra_host_obj_Float2D U,
                    chakra_host_obj_Float2D VT):
    """Call carma_svd_cpu
    
    chakra_host_obj_Float2D mat:
    chakra_host_obj_Float1D eigenvals:
    chakra_host_obj_Float2D U:
    chakra_host_obj_Float2D VT:
    """
    carma_svd_cpu[float](mat.c_h, eigenvals.c_h, VT.c_h, U.c_h)

def svd_Float(chakra_obj_Float2D mat,
                chakra_obj_Float1D eigenvals,
                chakra_obj_Float2D U,
                chakra_obj_Float2D VT):
    """Call carma_svd
    """

    carma_svd[float](mat.c_o, eigenvals.c_o, VT.c_o, U.c_o)


def getri_Float(chakra_obj_Float2D d_mat):
    return carma_getri[float](d_mat.c_o)


def getri_host_Float(chakra_host_obj_Float2D h_mat):
    cdef const long *dims=h_mat.c_h.getDims()
    return carma_getri_cpu[float](dims[1],h_mat.data_h)


def potri_Float(chakra_obj_Float2D d_mat ):
    return carma_potri[float](d_mat.c_o)


def potri_host_Float(chakra_host_obj_Float2D h_mat):
    cdef const long *dims=h_mat.c_h.getDims()
    return carma_potri_cpu[float](dims[1],h_mat.c_h.getData())


def syevd_Float(chakra_obj_Float2D d_A,
                np.ndarray[ndim=1,dtype=np.float32_t] eigenvals,
                chakra_obj_Float2D U=None,
                bool computeU=True):

    cdef long N=eigenvals.shape[0]
    
    if(U is None):
        if(computeU==True):
            carma_syevd[float]('V',N,d_A.c_o.getData(), <float*>eigenvals.data)
        else:
            carma_syevd[float]('N',N,d_A.c_o.getData(), <float*>eigenvals.data)
    else:
        U.c_o.copy(d_A.c_o,1,1)
        if(computeU==True):
            carma_syevd[float]('V',N,U.c_o.getData(), <float*>eigenvals.data)
        else:
            carma_syevd[float]('N',N,U.c_o.getData(), <float*>eigenvals.data)



def syevd_host_Float(chakra_host_obj_Float2D h_A,
                np.ndarray[ndim=1,dtype=np.float32_t] eigenvals,
                chakra_host_obj_Float2D U=None,
                bool computeU=True):
    
    cdef long N=eigenvals.shape[0]

    if(U is None):
        if(computeU==True):
            carma_syevd_cpu[float]('V',N, h_A.data_h, <float*>eigenvals.data)
        else:
            carma_syevd_cpu[float]('N',N, h_A.data_h, <float*>eigenvals.data)
    else:
        U.fill_from(h_A)
        if(computeU==True):
            carma_syevd_cpu[float]('V',N,U.data_h, <float*>eigenvals.data)
        else:
            carma_syevd_cpu[float]('N',N,U.data_h, <float*>eigenvals.data)




