/**
 * \class yoga_stream
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
#ifndef _YOGA_STREAM_H_
#define _YOGA_STREAM_H_

#include <yoga_utils.h>
#include <driver_types.h>
#include <vector>

using namespace std;

class yoga_streams {
protected:
	vector<cudaStream_t> streams;
	vector<cudaEvent_t> events;
	int eventflags;

public:
	yoga_streams();
	yoga_streams(unsigned int nbStreams);
	//yoga_stream(const yoga_stream& src_yoga_stream);
	~yoga_streams();

	cudaStream_t get_stream(int stream);
	cudaEvent_t get_event(int stream);
	cudaStream_t operator[](int idx) {
		return get_stream(idx);
	};

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

#endif // _YOGA_STREAM_H_
