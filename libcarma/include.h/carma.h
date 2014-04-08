/** 
 * \defgroup libcarma Carma
 *
 * \brief Carma is a library that provides GPU acceleration to Yorick
 *
 * \author $Author: Damien Gratadour & Arnaud Sevin $
 *
 */

#ifndef _CARMA_H_
#define _CARMA_H_

#include <carma_utils.h>
#include <carma_context.h>

//#define DEBUG 1
#define Y_SCOMPLEX 99

using namespace std;

extern "C" {
// TOOLS
void _carmaThreadExit();
void _carmaThreadSync();
carma_context *_getCurrentContext();

}
#endif // _CARMA_H_
