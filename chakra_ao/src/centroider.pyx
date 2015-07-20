cdef class Param_centroider:
    def set_type(self, bytes b):
        self.type_centro=b
    def set_type_fct(self, bytes b):
        self.type_fct=b
    def set_nwfs(self, long n):
        self.nwfs=n
    def set_nmax(self, long n):
        self.nmax=n
    def set_thresh(self, float t):
        self.thresh=t
    def set_width(self, float w):
        self.width=w
    def set_sizex(self, long s):
        self.sizex=s
    def set_sizey(self, long s):
        self.sizey=s
