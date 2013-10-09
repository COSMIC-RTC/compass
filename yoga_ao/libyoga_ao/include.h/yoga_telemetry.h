/**
 * \class yoga_telemetry
 *
 * \ingroup libyoga
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
#ifndef _YOGA_TELEMETRY_H_
#define _YOGA_TELEMETRY_H_

#include <yoga.h>
#include <yoga_obj.h>
#include <yoga_host_obj.h>
#include <yoga_streams.h>
#include <map>

typedef std::pair<std::string,int> type_telemetry_pair;

using namespace std;

class yoga_telemetry {
 protected:
  yoga_streams *streams;
  map< type_telemetry_pair, yoga_host_obj<float>* > objs;

 public:
  yoga_telemetry();
  yoga_telemetry(type_telemetry_pair obj,yoga_host_obj<float> *host_obj, unsigned int nbStreams);
  yoga_telemetry(std::string type_obj,int num_obj,yoga_host_obj<float> *host_obj, unsigned int nbStreams);
  //yoga_telemetry(const yoga_telemetry& src_yoga_telemetry);
  ~yoga_telemetry();

  //const yoga_host_obj<float>& operator[](int idx) const {return get_yoga_host_obj(idx);} ;

  int get_nbStreams();
  int add_stream();
  int add_stream(int nb);
  int del_stream();
  int del_stream(int nb);
  int del_all_streams();
  cudaStream_t get_cudaStream_t(int stream);

  int get_nbObjs();
  int add_obj(type_telemetry_pair obj, yoga_host_obj<float> *host_obj);
  int add_obj(std::string type_obj,int num_obj, yoga_host_obj<float> *host_obj);
  int del_obj(type_telemetry_pair obj);
  int del_obj(std::string type_obj,int num_obj);
  int del_all_objs();
  yoga_host_obj<float>* get_yoga_host_obj(type_telemetry_pair obj);
  yoga_host_obj<float>* get_yoga_host_obj(std::string type_obj,int num_obj);
  //yoga_host_obj<float>* operator[](int idx) {return get_yoga_host_obj(idx);} ;
  int cpy_obj(std::string type_obj,int num_obj, yoga_obj<float> *d_obj, cudaMemcpyKind flag);
  int wait_obj(std::string type_obj,int num_obj);

  int wait_stream(int stream);
  int wait_all_streams();

  /**< Memory transfer */
  int fill_from(std::string type_obj,int num_obj,float *data);
  int fill_into(std::string type_obj,int num_obj,float *data);

};

#endif // _YOGA_TELEMETRY_H_
