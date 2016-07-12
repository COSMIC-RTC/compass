/**
 * \class carma_stream
 *
 * \ingroup libcarma
 *
 * \brief this class
 *
 * \author $Author: dg, as $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2011/01/28$
 *
 */
#ifndef _CARMA_STREAM_H_
#define _CARMA_STREAM_H_

#include <carma_utils.h>
#include <driver_types.h>
#include <vector>

class carma_streams {
protected:
  std::vector<cudaStream_t> streams;
  std::vector<cudaEvent_t> events;
  int eventflags;

public:
  carma_streams();
  carma_streams(unsigned int nbStreams);
  //carma_stream(const carma_stream& src_carma_stream);
  ~carma_streams();

  cudaStream_t
  get_stream(int stream);
  cudaEvent_t
  get_event(int stream);
  cudaStream_t operator[](int idx) {
    return get_stream(idx);
  }

  int get_nbStreams();
  int add_stream();
  int add_stream(int nb);
  int del_stream();
  int del_stream(int nb);
  int del_all_streams();
  int wait_event(int stream);
  int wait_stream(int stream);
  int wait_all_streams();

};

#endif // _CARMA_STREAM_H_
