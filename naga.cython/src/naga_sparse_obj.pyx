

import numpy as np
cimport numpy as np
#np.import_array()

from scipy.sparse import csr_matrix

assert sizeof(int) == sizeof(np.int32_t)
assert sizeof(float) == sizeof(np.float32_t)
assert sizeof(double) == sizeof(np.float64_t)






##########################################################
#  P-Class `Float`
#########################################################
cdef class naga_sparse_obj_Float:

    def __cinit__(self,
                  naga_sparse_obj_Float obj=None,
                  np.ndarray[ndim=1, dtype=np.float32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Float constructor.

      carma_sparse_obj(carma_context *current_context)

      notimplemented
      carma_sparse_obj(carma_obj[float]* M)
      notimplemented
      carma_sparse_obj(carma_sparse_obj[float]* M)
      notimplemented
      #carma_sparse_obj(carma_context *current_context, carma_sparse_host_obj[float]* M)
      notimplemented
      carma_sparse_obj(carma_context *current_context, const long *dims, T_data * M, bool loadFromHost)
      notimplemented
      carma_sparse_obj(carma_context *current_context, const long *dims,
          T_data *values, int *colind, int *rowind, int nz, bool loadFromHost)

        input data must be a np.ndarray with ndim=1
        """

        cdef carma_context *ctxt=&carma_context.instance()

        if(obj is not None):
            #copy constructor
            #TODO
            ctxt.set_activeDevice(obj.getDevice(),1)
            self.c_sparse_obj=new carma_sparse_obj[float](obj.c_sparse_obj)
            return

        cdef int i

        if(data is None):
            if(dims is None):
                print("new carma_sparse_obj[float](ctxt)")
                self.c_sparse_obj=new carma_sparse_obj[float](ctxt)
                return

            if( dims.shape[0]!=1 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 1)")

            raise NotImplementedError
            #return

        #TODO
        cdef np.ndarray data_F=np.asfortranarray(data)
        raise NotImplementedError


    def __dealloc__(self):
        del self.c_sparse_obj



    cdef copy(self, carma_sparse_obj[float] *c_sparse_obj):
        del self.c_sparse_obj
        self.c_sparse_obj=new carma_sparse_obj[float](c_sparse_obj)



    def get_sparse(self):
        cdef int dims[2]
        dims[0]=self.c_sparse_obj.getDims(1)
        dims[1]=self.c_sparse_obj.getDims(2)
        cdef int nnz=self.c_sparse_obj.nz_elem

        cdef np.ndarray[ndim=1,dtype=np.int32_t]          rowInd=np.zeros(dims[0]+1,dtype=np.int32)
        cdef np.ndarray[ndim=1,dtype=np.int32_t]          colInd=np.zeros(nnz    ,dtype=np.int32)
        cdef np.ndarray[ndim=1,dtype=np.float32_t]    data=np.zeros(nnz,dtype=np.float32)


        self.c_sparse_obj.sparse_to_host(<int*>rowInd.data,<int*>colInd.data,<float*>data.data)
        return csr_matrix((data,colInd,rowInd),shape=(dims[0],dims[1]))




##########################################################
#  P-Class `Double`
#########################################################
cdef class naga_sparse_obj_Double:

    def __cinit__(self,
                  naga_sparse_obj_Double obj=None,
                  np.ndarray[ndim=1, dtype=np.float64_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Double constructor.

      carma_sparse_obj(carma_context *current_context)

      notimplemented
      carma_sparse_obj(carma_obj[double]* M)
      notimplemented
      carma_sparse_obj(carma_sparse_obj[double]* M)
      notimplemented
      #carma_sparse_obj(carma_context *current_context, carma_sparse_host_obj[double]* M)
      notimplemented
      carma_sparse_obj(carma_context *current_context, const long *dims, T_data * M, bool loadFromHost)
      notimplemented
      carma_sparse_obj(carma_context *current_context, const long *dims,
          T_data *values, int *colind, int *rowind, int nz, bool loadFromHost)

        input data must be a np.ndarray with ndim=1
        """

        cdef carma_context *ctxt=&carma_context.instance()

        if(obj is not None):
            #copy constructor
            #TODO
            ctxt.set_activeDevice(obj.getDevice(),1)
            self.c_sparse_obj=new carma_sparse_obj[double](obj.c_sparse_obj)
            return

        cdef int i

        if(data is None):
            if(dims is None):
                print("new carma_sparse_obj[double](ctxt)")
                self.c_sparse_obj=new carma_sparse_obj[double](ctxt)
                return

            if( dims.shape[0]!=1 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 1)")

            raise NotImplementedError
            #return

        #TODO
        cdef np.ndarray data_F=np.asfortranarray(data)
        raise NotImplementedError


    def __dealloc__(self):
        del self.c_sparse_obj



    cdef copy(self, carma_sparse_obj[double] *c_sparse_obj):
        del self.c_sparse_obj
        self.c_sparse_obj=new carma_sparse_obj[double](c_sparse_obj)



    def get_sparse(self):
        cdef int dims[2]
        dims[0]=self.c_sparse_obj.getDims(1)
        dims[1]=self.c_sparse_obj.getDims(2)
        cdef int nnz=self.c_sparse_obj.nz_elem

        cdef np.ndarray[ndim=1,dtype=np.int32_t]          rowInd=np.zeros(dims[0]+1,dtype=np.int32)
        cdef np.ndarray[ndim=1,dtype=np.int32_t]          colInd=np.zeros(nnz    ,dtype=np.int32)
        cdef np.ndarray[ndim=1,dtype=np.float64_t]    data=np.zeros(nnz,dtype=np.float64)


        self.c_sparse_obj.sparse_to_host(<int*>rowInd.data,<int*>colInd.data,<double*>data.data)
        return csr_matrix((data,colInd,rowInd),shape=(dims[0],dims[1]))


