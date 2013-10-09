/*
 * yoga_telemetry.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: sevin
 */

#include <yoga_telemetry.h>


yoga_telemetry::yoga_telemetry()
{
	streams = new yoga_streams();
}


yoga_telemetry::yoga_telemetry(type_telemetry_pair obj,yoga_host_obj<float> *host_obj, unsigned int nbStreams)
{
	add_obj(obj, host_obj);
	streams = new yoga_streams(nbStreams);
  //cutilSafeCall(cudaEventCreateWithFlags(&(this->start_event), cudaDeviceBlockingSync));
  //cutilSafeCall( cudaEventCreateWithFlags(&(this->stop_event), cudaDeviceBlockingSync));
  //cudaEventDefault
}

yoga_telemetry::yoga_telemetry(std::string type_obj,int num_obj,yoga_host_obj<float> *host_obj, unsigned int nbStreams)
{
	add_obj(type_obj, num_obj, host_obj);
	streams = new yoga_streams(nbStreams);
  //cutilSafeCall(cudaEventCreateWithFlags(&(this->start_event), cudaDeviceBlockingSync));
  //cutilSafeCall( cudaEventCreateWithFlags(&(this->stop_event), cudaDeviceBlockingSync));
  //cudaEventDefault
}

/*
yoga_telemetry::yoga_telemetry(const yoga_telemetry& src_yoga_telemetry)
{
	int nbStreams=src_yoga_telemetry.get_nbStreams();
	yoga_telemetry(nbStreams);
}
*/


yoga_telemetry::~yoga_telemetry()
{
	for( std::map< type_telemetry_pair, yoga_host_obj<float>* >::iterator it = objs.begin();
			objs.end()!= it ; it++) {
		delete it->second;
	}
	this->objs.clear();
	delete this->streams;

  //cudaEventDestroy(this->start_event);
  //cudaEventDestroy(this->stop_event);
}


int yoga_telemetry::get_nbObjs()
{
	return this->objs.size();
}


int yoga_telemetry::get_nbStreams()
{
	return this->streams->get_nbStreams();
}

int yoga_telemetry::add_stream()
{
	this->streams->add_stream();
	return this->streams->get_nbStreams();
}

int yoga_telemetry::add_stream(int nb)
{
	this->streams->add_stream(nb);
	return this->streams->get_nbStreams();
}

int yoga_telemetry::del_stream()
{
	this->streams->del_stream();
	return this->streams->get_nbStreams();
}

int yoga_telemetry::del_stream(int nb)
{
	this->streams->del_stream(nb);
	return this->streams->get_nbStreams();
}

int yoga_telemetry::del_all_streams()
{
	this->streams->del_all_streams();
	return this->streams->get_nbStreams();
}

int yoga_telemetry::add_obj(type_telemetry_pair obj,yoga_host_obj<float> *host_obj)
{
	this->objs[obj] = host_obj;
	return get_nbObjs();
}

int yoga_telemetry::add_obj(std::string type_obj,int num_obj,yoga_host_obj<float> *host_obj)
{
	this->objs[make_pair(type_obj, num_obj)] = host_obj;
	return get_nbObjs();
}

int yoga_telemetry::del_obj(type_telemetry_pair obj)
{
	delete this->objs[obj];
	this->objs.erase(obj);
	return get_nbObjs();
}

int yoga_telemetry::del_obj(std::string type_obj,int num_obj)
{
	delete this->objs[make_pair(type_obj, num_obj)];
	this->objs.erase(make_pair(type_obj, num_obj));
	return get_nbObjs();
}

yoga_host_obj<float>* yoga_telemetry::get_yoga_host_obj(type_telemetry_pair obj)
{
	return this->objs[obj];
}

yoga_host_obj<float>* yoga_telemetry::get_yoga_host_obj(std::string type_obj,int num_obj)
{
	return this->objs[make_pair(type_obj, num_obj)];
}

int yoga_telemetry::cpy_obj(std::string type_obj,int num_obj,yoga_obj<float> *d_obj, cudaMemcpyKind flag)
{
	this->objs[make_pair(type_obj, num_obj)]->cpy_obj(d_obj, flag);
	return EXIT_SUCCESS;
}

/**< Memory transfer */
int yoga_telemetry::fill_from(std::string type_obj,int num_obj,float *data)
{
	this->objs[make_pair(type_obj, num_obj)]->fill_from(data);
	return EXIT_SUCCESS;
}

int yoga_telemetry::fill_into(std::string type_obj,int num_obj,float *data)
{
	this->objs[make_pair(type_obj, num_obj)]->fill_into(data);
	return EXIT_SUCCESS;
}

int yoga_telemetry::wait_obj(std::string type_obj,int num_obj){
	this->objs[make_pair(type_obj, num_obj)]->wait_all_streams();
	return EXIT_SUCCESS;
}

cudaStream_t yoga_telemetry::get_cudaStream_t(int stream)
{
	return this->streams->get_stream(stream);
}


int yoga_telemetry::wait_stream(int stream)
{
	this->streams->wait_stream(stream);
	return EXIT_SUCCESS;
}


int yoga_telemetry::wait_all_streams()
{
	this->streams->wait_all_streams();
	return EXIT_SUCCESS;
}
