#include <carma_streams.h>
#include <carma_utils.h>

carma_streams::carma_streams()
{
	carma_streams(0);
}

carma_streams::carma_streams(unsigned int nbStreams)
{
  //this->streams  = new vector<cudaStream_t>();
  for(unsigned int i = 0; i < nbStreams; i++) {
	  add_stream();
  }

  //cutilSafeCall(cudaEventCreateWithFlags(&(this->start_event), cudaDeviceBlockingSync));
  //cutilSafeCall( cudaEventCreateWithFlags(&(this->stop_event), cudaDeviceBlockingSync));
  //cudaEventDefault
}

/*
carma_stream::carma_stream(const carma_stream& src_carma_stream)
{
	int nbStreams=src_carma_stream.get_nbStreams();
	carma_stream(nbStreams);
}
*/

carma_streams::~carma_streams()
{
  del_all_streams();
  this->streams.clear();

  //cudaEventDestroy(this->start_event);
  //cudaEventDestroy(this->stop_event);
}

int carma_streams::get_nbStreams()
{
	return this->streams.size();
}

int carma_streams::add_stream()
{
	  cudaStream_t stream_tmp;
	  cutilSafeCall( cudaStreamCreate(&stream_tmp) );
	  this->streams.push_back(stream_tmp);

	  cudaEvent_t event_tmp;
	  cutilSafeCall( cudaEventCreate(&event_tmp) );
	  this->events.push_back(event_tmp);

#ifdef DEBUG
	  printf("CARMA Stream created @ %8.8lX\n", (unsigned long)stream_tmp);
#endif

	  return get_nbStreams();
}

int carma_streams::add_stream(int nb)
{
	for(int stream=0; stream<nb; stream++)
		add_stream();
	return get_nbStreams();
}

int carma_streams::del_stream()
{
	if(streams.empty()) return 0;

	cudaStream_t stream_tmp = this->streams[this->streams.size()-1];
	cutilSafeCall( cudaStreamDestroy(stream_tmp) );
	this->streams.pop_back();

	cudaEvent_t event_tmp = this->events[this->events.size()-1];
	cutilSafeCall( cudaEventDestroy(event_tmp) );
	this->events.pop_back();

#ifdef DEBUG
	printf("CARMA Stream deleted @ %8.8lX\n", (unsigned long)stream_tmp);
#endif
	return get_nbStreams();
}

int carma_streams::del_stream(int nb)
{
	for(int stream=0; stream<nb && !streams.empty(); stream++)
		del_stream();
	return get_nbStreams();
}

int carma_streams::del_all_streams()
{
	while(!streams.empty())
		del_stream();
	return get_nbStreams();
}

cudaStream_t carma_streams::get_stream(int stream)
{
	return this->streams[stream];
}

cudaEvent_t carma_streams::get_event(int stream)
{
	return this->events[stream];
}

int carma_streams::wait_event(int stream)
{
	  cutilSafeCall( cudaEventSynchronize(this->events[stream]) );
	  return EXIT_SUCCESS;
}

int carma_streams::wait_stream(int stream)
{
	  cutilSafeCall( cudaStreamSynchronize(this->streams[stream]) );
	  return EXIT_SUCCESS;
}

int carma_streams::wait_all_streams()
{
	for(unsigned int stream=0; stream<streams.size(); stream++)
	  cutilSafeCall( cudaStreamSynchronize(this->streams[stream]) );
	return EXIT_SUCCESS;
}

