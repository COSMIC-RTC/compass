/*! @file yoga_api.h

 @brief Definition of functions callable by the Yorick interpreter

 Based on the standard API for interfacing Yorick packages to the interpreter
 Binding is realized through the definition of the corresponding "extern" wrappers
 in the yoga.i file (see yorick directory). Functions callable by the interpreter
 start with a Y_

 Several persistent objects are also defined :
 - yContext
 - yoga_yObj
 - yoga_host_yObj

 Copyright (c) 2011, Damien Gratadour & Arnaud Sevin

 This file is part of YoGA, the Yorick with GPU Acceleration plugin

 This program is free software; you can redistribute it and/or  modify it
 under the terms of the GNU General Public License  as  published  by the
 Free Software Foundation; either version 3 of the License,  or  (at your
 option) any later version.
 
 This program is distributed in the hope  that  it  will  be  useful, but
 WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _YOGA_API_H_
#define _YOGA_API_H_

#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_sparse_host_obj.h>
#include <carma_sparse_obj.h>
#include <carma_cublas.h>
#include <carma_cusparse.h>
#include <carma_utils.h>
#include <carma_ipcs.h>
#include <carma_timer.h>
#include <carma.h>

#include <yapi.h>
#define Y_SCOMPLEX 99

extern "C" {

/*    _                        _       __
 *   | |_ _   _ _ __   ___  __| | ___ / _|
 *   | __| | | | '_ \ / _ \/ _` |/ _ \ |_
 *   | |_| |_| | |_) |  __/ (_| |  __/  _|
 *    \__|\__, | .__/ \___|\__,_|\___|_|
 *        |___/|_|
 */
typedef struct context_struct {
  /** 
   * @typedef Yorick API context structure
   */
  void *carma_context; /**< pointer to a carma_context object */
} context_struct;


typedef struct timer_struct {
  /**
   * @typedef Yorick API context structure
   */
  void *carma_timer; /**< pointer to a carma_context object */
} timer_struct;



typedef struct ipcs_struct {
  /**
   * @typedef Yorick API ipcs structure
   */
  void *carma_ipcs; /**< pointer to a carma_ipcs object */
} ipcs_struct;



typedef struct yObj_struct {
  /** 
   * @typedef Yorick API yoga_object structure
   */
  void *carma_object; /**< pointer to a yoga_object */
  int type; /**< type of data in the yoga_object (Yorick API types) */
  int device; /**< device number on which it resides */
  unsigned char isRef; /**< reference counting */
} yObj_struct;

typedef struct yHostObj_struct {
  /** 
   * @typedef Yorick API yoga_host_object structure
   */
  void *carma_host_object; /**< pointer to a yoga_host_object */
  int type; /**< type of data in the yoga_host_object (Yorick API types) */
} yHostObj_struct;

typedef struct ySparseObj_struct {
  /**
   * @typedef Yorick API yoga_object structure
   */
  void *carma_sparse_object; /**< pointer to a yoga_sparse_object */
  int type; /**< type of data in the yoga_sparse_object (Yorick API types) */
  int device; /**< device number on which it resides */
  unsigned char isRef; /**< reference counting */
} ySparseObj_struct;

typedef struct ySparseHostObj_struct {
  /**
   * @typedef Yorick API yoga_host_object structure
   */
  void *carma_sparse_host_object; /**< pointer to a yoga_sparse_host_object */
  int type; /**< type of data in the yoga_sparse_host_object (Yorick API types) */
} ySparseHostObj_struct;

/*
 */

void
timer_print(void *obj);
void
timer_free(void *obj);

timer_struct*
yoga_getTimer(int argc, int pos);

void
Y_yoga_timer(int argc);

void
Y_yoga_timer_start(int argc);

void
Y_yoga_timer_stop(int argc);

void
Y_yoga_timer_reset(int argc);

/*
 *  ___ ____   ____ ____
 * |_ _|  _ \ / ___/ ___|
 *  | || |_) | |   \___ \
 *  | ||  __/| |___ ___) |
 * |___|_|    \____|____/
*/

void
ipcs_print(void *obj);
void
ipcs_free(void *obj);

ipcs_struct*
yoga_getIPCs(int argc, int pos);

void
Y_yoga_ipcs(int argc);

//register an event, id must be a non nul argument
//int Y_yoga_register_cuevent(unsigned int id, CUevent event);

//get a memory handle
//int Y_yoga_get_memHandle(unsigned int id, CUipcMemHandle *phandle);
//get a event handle
//int Y_yoga_get_eventHandle(unsigned int id, CUipcEventHandle *phandle);

//free a memory handle shared mem space
void Y_yoga_free_memHandle(int argc);
//free a event handle shared event space
//void Y_yoga_free_eventHandle(unsigned int id);

/*
  Transfer via CPU memory methods
*/
//allocation of the shm for memory tranfers
//int Y_yoga_alloc_memtransfer_shm(unsigned int id, void *shm, size_t bsize);
//get tranfer shm
//int Y_yoga_get_memtransfer_shm(unsigned int id, void *shm);
//free transfer shm ref by id
//void Y_yoga_free_memtransfer_shms(unsigned int id);


/*
  Barrier methods
*/
void Y_yoga_init_barrier(int argc);

void Y_yoga_wait_barrier(int argc);

void Y_yoga_free_barrier(int argc);


/*                  _            _
 *   ___ ___  _ __ | |_ _____  _| |_
 *  / __/ _ \| '_ \| __/ _ \ \/ / __|
 * | (_| (_) | | | | ||  __/>  <| |_
 *  \___\___/|_| |_|\__\___/_/\_\\__|
 */

void
context_print(void *obj);
void
context_free(void *obj);

carma_context*
_getCurrentContext();

context_struct*
yoga_getContext(int argc, int pos);

void
Y_yoga_context(int argc);
void
Y_context_getactivedevice(int argc);
void
Y_context_get_device_info(int argc);
void
Y_context_get_device_meminfo(int argc);
void
Y_context_get_maxGflopsDeviceId(int argc);
void
Y_activeDevice(int argc);
void
_yoga_setDevice(int mydevice);
int
_yoga_getDevice();
int
_yoga_getnDevice();
void
_yogaThreadExit();
void
_yogaThreadSync();
void
_yoga_init();

/*                                    _     _           _
 *   _   _  ___   __ _  __ _     ___ | |__ (_) ___  ___| |_
 *  | | | |/ _ \ / _` |/ _` |   / _ \| '_ \| |/ _ \/ __| __|
 *  | |_| | (_) | (_| | (_| |  | (_) | |_) | |  __/ (__| |_
 *   \__, |\___/ \__, |\__,_|___\___/|_.__// |\___|\___|\__|
 *   |___/       |___/     |_____|       |__/
 */

void
yObj_print(void *obj);
void
yObj_eval(void *obj, int n);
void
yObj_free(void *obj);

yObj_struct*
yoga_getyObj(int argc, int pos);

void
Y_yoga_obj(int argc);
void
Y_yoga_getp(int argc);
void
Y_yoga_setv(int argc);
void
Y_yoga_setm(int argc);
void
Y_yoga_host2device(int argc);
void
Y_yoga_device2host(int argc);
void
Y_yoga_imin(int argc);
void
Y_yoga_imax(int argc);
void
Y_yoga_asum(int argc);
void
Y_yoga_nrm2(int argc);
void
Y_yoga_scale(int argc);
void
Y_yoga_swap(int argc);
void
Y_yoga_axpy(int argc);
void
Y_yoga_dot(int argc);
void
Y_yoga_mv(int argc);
void
Y_yoga_rank1(int argc);
void
Y_yoga_mm(int argc);
void
Y_yoga_transpose(int argc);
void
Y_yoga_random(int argc);
void
Y_yoga_random_n(int argc);
void
Y_yoga_fft(int argc);
void
Y_yoga_max(int argc);
void
Y_yoga_min(int argc);
void
Y_yoga_mult(int argc);
void
Y_yoga_add(int argc);
void
Y_yoga_sort(int argc);
void
Y_yoga_compact(int argc);
void
Y_yoga_fftconv(int argc);
void
Y_yoga_fftconv_init(int argc);

/*                               _               _
 *   _   _  ___   __ _  __ _    | |__   ___  ___| |_
 *  | | | |/ _ \ / _` |/ _` |   | '_ \ / _ \/ __| __|
 *  | |_| | (_) | (_| | (_| |   | | | | (_) \__ \ |_
 *   \__, |\___/ \__, |\__,_|___|_| |_|\___/|___/\__|
 *   |___/       |___/     |_____|
 */

void
yHostObj_free(void *obj);
void
yHostObj_print(void *obj);
void
yHostObj_eval(void *obj, int n);

void
Y_yoga_host_obj(int argc);
void
Y_yoga_host_getp(int argc);

/*                                                                     _     _
 *  _   _  ___   __ _  __ _     ___ _ __   __ _ _ __ ___  ___     ___ | |__ (_)
 * | | | |/ _ \ / _` |/ _` |   / __| '_ \ / _` | '__/ __|/ _ \   / _ \| '_ \| |
 * | |_| | (_) | (_| | (_| |   \__ \ |_) | (_| | |  \__ \  __/  | (_) | |_) | |
 *  \__, |\___/ \__, |\__,_|___|___/ .__/ \__,_|_|  |___/\___|___\___/|_.__// |
 *  |___/       |___/     |_____|  |_|                      |_____|       |__/
 */

void
ySparseObj_print(void *obj);
void
ySparseObj_eval(void *obj, int n);
void
ySparseObj_free(void *obj);

void
Y_yoga_sparse_obj(int argc);

/*                                                                _               _            _     _
 *  _   _  ___   __ _  __ _     ___ _ __   __ _ _ __ ___  ___    | |__   ___  ___| |_     ___ | |__ (_)
 * | | | |/ _ \ / _` |/ _` |   / __| '_ \ / _` | '__/ __|/ _ \   | '_ \ / _ \/ __| __|   / _ \| '_ \| |
 * | |_| | (_) | (_| | (_| |   \__ \ |_) | (_| | |  \__ \  __/   | | | | (_) \__ \ |_   | (_) | |_) | |
 *  \__, |\___/ \__, |\__,_|___|___/ .__/ \__,_|_|  |___/\___|___|_| |_|\___/|___/\__|___\___/|_.__// |
 *  |___/       |___/     |_____|  |_|                      |_____|                 |_____|       |__/
 */

void
ySparseHostObj_print(void *obj);
void
ySparseHostObj_eval(void *obj, int n);
void
ySparseHostObj_free(void *obj);

void
Y_yoga_sparse_host_obj(int argc);

/*    _ __ ___   __ _  __ _ _ __ ___   __ _
 *   | '_ ` _ \ / _` |/ _` | '_ ` _ \ / _` |
 *   | | | | | | (_| | (_| | | | | | | (_| |
 *   |_| |_| |_|\__,_|\__, |_| |_| |_|\__,_|
 *                    |___/
 */

void
Y_yoga_svd(int argc);
void
Y_yoga_svd_host(int argc);

/*             _
 *   ___ _   _| | __ _
 *  / __| | | | |/ _` |
 * | (__| |_| | | (_| |
 *  \___|\__,_|_|\__,_|
 */

void
Y_yoga_cula_svd(int argc);
void
Y_yoga_cula_svd_host(int argc);

/*    __ _ _ __ _ __ __ _ _   _ ___
 *   / _` | '__| '__/ _` | | | / __|
 *  | (_| | |  | | | (_| | |_| \__ \
 *   \__,_|_|  |_|  \__,_|\__, |___/
 *                        |___/
 */

void
Y_yoga_getarray(int argc);
void
Y_yoga_fillarray(int argc);
void
Y_yoga_getvalue(int argc);
void
Y_yoga_plus(int argc);
void
Y_yoga_plusai(int argc);
void
Y_yoga_test(int argc);

}

#endif // _YOGA_API_H_
