#include <yoga_streams.h>
#include <yoga_utils.h>

yoga_streams::yoga_streams()
{
	yoga_streams(0);
}

yoga_streams::yoga_streams(unsigned int nbStreams)
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
yoga_stream::yoga_stream(const yoga_stream& src_yoga_stream)
{
	int nbStreams=src_yoga_stream.get_nbStreams();
	yoga_stream(nbStreams);
}
*/

yoga_streams::~yoga_streams()
{
  del_all_streams();
  this->streams.clear();

  //cudaEventDestroy(this->start_event);
  //cudaEventDestroy(this->stop_event);
}

int yoga_streams::get_nbStreams()
{
	return this->streams.size();
}

int yoga_streams::add_stream()
{
	  cudaStream_t stream_tmp;
	  cutilSafeCall( cudaStreamCreate(&stream_tmp) );
	  this->streams.push_back(stream_tmp);

	  cudaEvent_t event_tmp;
	  cutilSafeCall( cudaEventCreate(&event_tmp) );
	  this->events.push_back(event_tmp);

#ifdef DEBUG
	  printf("YOGA Stream created @ %8.8lX\n", (unsigned long)stream_tmp);
#endif

	  return get_nbStreams();
}

int yoga_streams::add_stream(int nb)
{
	for(int stream=0; stream<nb; stream++)
		add_stream();
	return get_nbStreams();
}

int yoga_streams::del_stream()
{
	if(streams.empty()) return 0;

	cudaStream_t stream_tmp = this->streams[this->streams.size()-1];
	cutilSafeCall( cudaStreamDestroy(stream_tmp) );
	this->streams.pop_back();

	cudaEvent_t event_tmp = this->events[this->events.size()-1];
	cutilSafeCall( cudaEventDestroy(event_tmp) );
	this->events.pop_back();

#ifdef DEBUG
	printf("YOGA Stream deleted @ %8.8lX\n", (unsigned long)stream_tmp);
#endif
	return get_nbStreams();
}

int yoga_streams::del_stream(int nb)
{
	for(int stream=0; stream<nb && !streams.empty(); stream++)
		del_stream();
	return get_nbStreams();
}

int yoga_streams::del_all_streams()
{
	while(!streams.empty())
		del_stream();
	return get_nbStreams();
}

cudaStream_t yoga_streams::get_stream(int stream)
{
	return this->streams[stream];
}

cudaEvent_t yoga_streams::get_event(int stream)
{
	return this->events[stream];
}

int yoga_streams::wait_event(int stream)
{
	  cutilSafeCall( cudaEventSynchronize(this->events[stream]) );
	  return EXIT_SUCCESS;
}

int yoga_streams::wait_stream(int stream)
{
	  cutilSafeCall( cudaStreamSynchronize(this->streams[stream]) );
	  return EXIT_SUCCESS;
}

int yoga_streams::wait_all_streams()
{
	for(unsigned int stream=0; stream<streams.size(); stream++)
	  cutilSafeCall( cudaStreamSynchronize(this->streams[stream]) );
	return EXIT_SUCCESS;
}

