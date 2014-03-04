/**
 * \class sutra_telemetry
 *
 * \ingroup libsutra
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
#ifndef _SUTRA_TELEMETRY_H_
#define _SUTRA_TELEMETRY_H_

#include <carma.h>
#include <carma_obj.h>
#include <carma_host_obj.h>
#include <carma_streams.h>
#include <map>

typedef std::pair<std::string, int> type_telemetry_pair;

using namespace std;

class sutra_telemetry {
protected:
	carma_streams *streams;
	map<type_telemetry_pair, carma_host_obj<float>*> objs;

public:
	sutra_telemetry();
	sutra_telemetry(type_telemetry_pair obj, carma_host_obj<float> *host_obj,
			unsigned int nbStreams);
	sutra_telemetry(std::string type_obj, int num_obj,
			carma_host_obj<float> *host_obj, unsigned int nbStreams);
	//sutra_telemetry(const sutra_telemetry& src_sutra_telemetry);
	~sutra_telemetry();

	//const carma_host_obj<float>& operator[](int idx) const {return get_sutra_host_obj(idx);} ;

	int get_nbStreams();
	int add_stream();
	int add_stream(int nb);
	int del_stream();
	int del_stream(int nb);
	int del_all_streams();
	cudaStream_t get_cudaStream_t(int stream);

	int get_nbObjs();
	int add_obj(type_telemetry_pair obj, carma_host_obj<float> *host_obj);
	int add_obj(std::string type_obj, int num_obj,
			carma_host_obj<float> *host_obj);
	int del_obj(type_telemetry_pair obj);
	int del_obj(std::string type_obj, int num_obj);
	int del_all_objs();
	carma_host_obj<float>* get_carma_host_obj(type_telemetry_pair obj);
	carma_host_obj<float>* get_carma_host_obj(std::string type_obj,
			int num_obj);
	//carma_host_obj<float>* operator[](int idx) {return get_carma_host_obj(idx);} ;
	int cpy_obj(std::string type_obj, int num_obj, carma_obj<float> *d_obj,
			cudaMemcpyKind flag);
	int wait_obj(std::string type_obj, int num_obj);

	int wait_stream(int stream);
	int wait_all_streams();

	/**< Memory transfer */
	int fill_from(std::string type_obj, int num_obj, float *data);
	int fill_into(std::string type_obj, int num_obj, float *data);

};

#endif // _SUTRA_TELEMETRY_H_
