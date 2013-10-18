/** 
 * \defgroup libyoga Yoga
 *
 * \brief Yoga is a library that provides GPU acceleration to Yorick
 *
 * \author $Author: Damien Gratadour & Arnaud Sevin $
 *
 */

#ifndef _YOGA_H_
#define _YOGA_H_

#include <yoga_utils.h>
#include <yoga_context.h>

//#define DEBUG 1
#define Y_SCOMPLEX 99

using namespace std;

extern "C" {
  // TOOLS
  void _yogaThreadExit();
  void _yogaThreadSync();
  yoga_context *_getCurrentContext();

}
#endif // _YOGA_H_
