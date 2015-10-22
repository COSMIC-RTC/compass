from libcpp.vector cimport vector

cdef extern from "carma_streams.h":
    cdef cppclass carma_streams:
        int eventflags
        carma_streams()
        carma_streams(unsigned int nbStreams)
        #~carma_streams()

        cudaStream_t get_stream(int stream)
        cudaEvent_t get_event(int stream)
        cudaStream_t operator[](int idx)

        int get_nbStreams()
        int add_stream()
        int add_stream(int nb)
        int del_stream()
        int del_stream(int nb)
        int del_all_streams()
        int wait_event(int stream)
        int wait_stream(int stream)
        int wait_all_streams()
