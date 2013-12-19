/*! @file yoga.cpp

 @brief Code for functions callable by the Yorick interpreter

 see yoga_api.h for a complete description

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

#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_cublas.h>
#include <carma.h>
#include <sstream>
#include <iomanip>
#include <cula.hpp>
#include "yoga_api.h"
#include <cudaProfiler.h>

static y_userobj_t yContext =
  {
    /** 
     * @typedef Yorick API context userobj
     */
    const_cast<char*>("yContext Object"), &context_free, &context_print, 0, 0, 0
  };

static y_userobj_t yObj =
  {
    /** 
     * @typedef Yorick API object_yObj userobj
     */
    const_cast<char*>("Carma Object"), &yObj_free, &yObj_print, &yObj_eval, 0, 0
  };

static y_userobj_t yHostObj =
  {
    /** 
     * @typedef Yorick API host_obj userobj
     */
    const_cast<char*>("Carma Host Object"), &yHostObj_free, &yHostObj_print, &yHostObj_eval, 0, 0
  };

extern "C" {

/*
 %                  _            _
 *   ___ ___  _ __ | |_ _____  _| |_
 *  / __/ _ \| '_ \| __/ _ \ \/ / __|
 * | (_| (_) | | | | ||  __/>  <| |_
 *  \___\___/|_| |_|\__\___/_/\_\\__|
 */

void context_free(void *obj) {
	/** @brief context_struct destructor.
	 *  @param obj : context_struct to freed
	 */
	context_struct *handler = (context_struct *) obj;
	try {
		carma_context *context_obj_handler =
				(carma_context *) (handler->carma_context);
		delete context_obj_handler;
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	}
}

void context_print(void *obj)
/** @brief context_struct printer.
 *  @param[in] obj : context_struct to print
 */
{
	context_struct *handler = (context_struct *) obj;
	carma_context *context_handler = (carma_context *) (handler->carma_context);
	unsigned int activeDevice = context_handler->get_activeDevice();
	unsigned int nDevice = context_handler->get_ndevice();

	cout << "CArMA Context : " << endl;
	cout << nDevice << " CUDA device(s)" << endl;
	cout << "Dev Id" << " | " << "multiproc" << " | " << "cores" << " | "
			<< "compute" << " | " << "Perf (GFlops)" << endl;
	for (size_t idx = 0; idx < nDevice; idx++) {
		cout << ((idx == activeDevice) ? "<U>" : "   ") << setw(3) << idx
				<< " | " << setw(9)
				<< context_handler->get_device(idx)->get_properties().multiProcessorCount
				<< " | " << setw(5)
				<< context_handler->get_device(idx)->get_sm_per_multiproc()
						* context_handler->get_device(idx)->get_properties().multiProcessorCount
				<< " | " << setw(5)
				<< context_handler->get_device(idx)->get_properties().major
				<< "."
				<< context_handler->get_device(idx)->get_properties().minor
				<< " | "
				<< context_handler->get_device(idx)->get_compute_perf() / 1.e6
				<< endl;
	}
}

void Y_yoga_context(int argc)
/** @brief context_struct creator.
 *  @param[in] argc : command line argument (current context expected)
 */
{
	try {
		context_struct *handle = (context_struct *) ypush_obj(&yContext,
				sizeof(context_struct));
		handle->carma_context = new carma_context();
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	} catch (...) {
		stringstream buf;
		buf << "unknown error with carma_context construction in " << __FILE__
				<< "@" << __LINE__ << endl;
		y_error(buf.str().c_str());
	}
}

carma_context*
_getCurrentContext()
/** @brief simple routine to retrieve current context
 *  @return the current carma_context
 */
{
	ypush_global(yfind_global("current_context\0", 0));
	context_struct *handle = (context_struct *) yget_obj(0, &yContext);
	yarg_drop(1);
	return (carma_context *) handle->carma_context;
}

void Y_context_getactivedevice(int argc)
/** @brief simple routine to retrieve current active device
 *  @param[in] argc : command line argument (current context expected)
 *  the current active device id is pushed on the stack
 */
{
	try {
		context_struct *handle = (context_struct *) yget_obj(argc - 1,
				&yContext);
		carma_context *context_handle = (carma_context *) handle->carma_context;
		int activeDevice = context_handle->get_activeDevice();
		ypush_int(activeDevice);
	} catch (...) {
		stringstream buf;
		buf << "unknown error with carma_context in " << __FILE__ << "@"
				<< __LINE__ << endl;
		y_error(buf.str().c_str());
	}
}

void Y_context_get_maxGflopsDeviceId(int argc)
/** @brief simple routine to retrieve current active device
 *  @param[in] argc : command line argument (current context expected)
 *  the current active device id is pushed on the stack
 */
{
	try {
		context_struct *handle = (context_struct *) yget_obj(argc - 1,
				&yContext);
		carma_context *context_handle = (carma_context *) handle->carma_context;
		int device = context_handle->get_maxGflopsDeviceId();
		ypush_int(device);
	} catch (...) {
		stringstream buf;
		buf << "unknown error with carma_context in " << __FILE__ << "@"
				<< __LINE__ << endl;
		y_error(buf.str().c_str());
	}
}

void Y_activeDevice(int argc)
/** @brief simple routine to get / set a device as active
 *  callable as a subroutine or as a function
 *  if called as a subroutine, parameter gives the device id
 *  if called as a function, push the active device id on the stack
 *  @param[in] mydevice : the device id
 */
{
	carma_context *context_handle = _getCurrentContext();
	if (yarg_subroutine()) {
		int odevice = ygets_i(argc-1);
		context_handle->set_activeDevice(odevice, 0);
	} else {
		ypush_int(context_handle->get_activeDevice());
	}
}

void _yoga_setDevice(int mydevice)
/** @brief simple routine to set a device as active
 *  @param[in] mydevice : the device id
 */
{
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(mydevice);
}

int _yoga_getDevice()
/** @brief simple wrapper to get the active device
 */
{
	carma_context *context_handle = _getCurrentContext();
	return context_handle->get_activeDevice();
}

int _yoga_getnDevice()
/** @brief simple wrapper to get the number of usable devices
 */
{
	carma_context *context_handle = _getCurrentContext();
	return context_handle->get_ndevice();
}

void _yogaThreadExit()
/*! \brief simple wrapper for general ThreadExist
 */
{
	//cerr << "Shutting down : " ;
	culaShutdown();
	cutilSafeCall( cudaThreadExit());
	//cerr << "OK " << __LINE__ << endl;
}

void _yogaThreadSync()
/*! \brief simple wrapper for general threads synchronization
 */
{
	cutilSafeThreadSync();
}

void _yoga_init()
/*! \brief simple routine for general init at YoGA launch
 */
{
	culaStatus status = culaInitialize();
	if (status) {
		char buf[256];
		culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
		printf("%s\n", buf);
	}
	ycall_on_quit(_yogaThreadExit);
}

void _yoga_start_profiler() {
	  cuProfilerStart();
}
void _yoga_stop_profiler() {
	  cuProfilerStop();
}

/*
 *                                   _     _           _
 *  _   _  ___   __ _  __ _     ___ | |__ (_) ___  ___| |_
 * | | | |/ _ \ / _` |/ _` |   / _ \| '_ \| |/ _ \/ __| __|
 * | |_| | (_) | (_| | (_| |  | (_) | |_) | |  __/ (__| |_
 *  \__, |\___/ \__, |\__,_|___\___/|_.__// |\___|\___|\__|
 *  |___/       |___/     |_____|       |__/
 */

void yObj_print(void *obj)
/** @brief yObj_struct printer.
 *  @param[in] obj : yObj_struct to print
 */
{
	ostringstream mystr;
	yObj_struct *handler = (yObj_struct *) obj;
	try {
		if (handler->type == Y_FLOAT) {
			carma_obj<float> *carma_obj_handler =
					(carma_obj<float> *) (handler->carma_object);
			mystr << "Carma Object : " << endl;
			mystr << "type : float" << endl;
			long *dims = carma_obj_handler->getDims();
			mystr << "Dimensions : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_DOUBLE) {
			carma_obj<double> *carma_obj_handler =
					(carma_obj<double> *) (handler->carma_object);
			mystr << "Carma Object : " << endl;
			mystr << "type : double" << endl;
			long *dims = carma_obj_handler->getDims();
			mystr << "Dimensions : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_INT) {
			carma_obj<int> *carma_obj_handler =
					(carma_obj<int> *) (handler->carma_object);
			mystr << "Carma Object : " << endl;
			mystr << "type : int" << endl;
			long *dims = carma_obj_handler->getDims();
			mystr << "Dimensions : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_SHORT) {
			carma_obj<unsigned int> *carma_obj_handler =
					(carma_obj<unsigned int> *) (handler->carma_object);
			mystr << "Carma Object : " << endl;
			mystr << "type : unsigned int" << endl;
			long *dims = carma_obj_handler->getDims();
			mystr << "Dimensions : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_COMPLEX) {
			carma_obj<cuDoubleComplex> *carma_obj_handler = (carma_obj<
					cuDoubleComplex> *) (handler->carma_object);
			mystr << "Carma Object : " << endl;
			mystr << "type : double complex" << endl;
			long *dims = carma_obj_handler->getDims();
			mystr << "Dimensions : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_SCOMPLEX) {
			carma_obj<cuDoubleComplex> *carma_obj_handler = (carma_obj<
					cuDoubleComplex> *) (handler->carma_object);
			mystr << "Carma Object : " << endl;
			mystr << "type : simple complex" << endl;
			long *dims = carma_obj_handler->getDims();
			mystr << "Dimensions : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else
			throw "Type unknown";
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	}
	y_print(mystr.str().c_str(), 1);
}

void yObj_eval(void *obj, int n)
/** @brief yObj_struct evaluator.
 *  @param[in] obj : yObj_struct to evaluate
 */
{
	yObj_struct *handler = (yObj_struct *) obj;

	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDeviceForCpy(handler->device);

	try {
		if (handler->type == Y_FLOAT) {
			caObjS *carma_obj_handler = (caObjS *) (handler->carma_object);
			float *data = ypush_f(carma_obj_handler->getDims());
			carma_obj_handler->device2host(data);
		} else if (handler->type == Y_DOUBLE) {
			caObjD *carma_obj_handler = (caObjD *) (handler->carma_object);
			double *data = ypush_d(carma_obj_handler->getDims());
			carma_obj_handler->device2host(data);
		} else if (handler->type == Y_INT) {
			caObjI *carma_obj_handler = (caObjI *) (handler->carma_object);
			int *data = ypush_i(carma_obj_handler->getDims());
			carma_obj_handler->device2host(data);
		} else if (handler->type == Y_SHORT) {
			caObjUI *carma_obj_handler = (caObjUI *) (handler->carma_object);
			unsigned int *data = (unsigned int *) ypush_i(
					carma_obj_handler->getDims());
			carma_obj_handler->device2host(data);
		} else if (handler->type == Y_COMPLEX) {
			caObjZ *carma_obj_handler = (caObjZ *) (handler->carma_object);
			cuDoubleComplex *data = (cuDoubleComplex *) ypush_z(
					carma_obj_handler->getDims());
			carma_obj_handler->device2host(data);
		} else if (handler->type == Y_SCOMPLEX) {
			caObjC *carma_obj_handler = (caObjC *) (handler->carma_object);
			/* scomplex -> float2 */
			long int *ndims_obj = carma_obj_handler->getDims();
			long int *ndims_data = new long[ndims_obj[0] + 2];
			ndims_data[0] = ndims_obj[0] + 1;
			ndims_data[1] = 2;
			memcpy(&ndims_data[2], &ndims_obj[1], sizeof(long) * ndims_obj[0]);
			float *h_data = (float *) ypush_f(ndims_data);
			carma_obj_handler->device2host((cuFloatComplex *) h_data);
			delete ndims_data;
		} else
			throw "Type unknown";
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	}
}

void yObj_free(void *obj)
/** @brief yObj_struct destructor.
 *  @param[in] obj : yObj_struct to free
 */
{
	yObj_struct *handler = (yObj_struct *) obj;
	if (handler->isRef)
		return;
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(handler->device);
	try {
		if (handler->type == Y_FLOAT) {
			caObjS *carma_obj_handler = (caObjS *) (handler->carma_object);
			delete carma_obj_handler;
		} else if (handler->type == Y_DOUBLE) {
			caObjD *carma_obj_handler = (caObjD *) (handler->carma_object);
			delete carma_obj_handler;
		} else if (handler->type == Y_INT) {
			caObjI *carma_obj_handler = (caObjI *) (handler->carma_object);
			delete carma_obj_handler;
		} else if (handler->type == Y_SHORT) {
			caObjUI *carma_obj_handler = (caObjUI *) (handler->carma_object);
			delete carma_obj_handler;
		} else if (handler->type == Y_COMPLEX) {
			caObjZ *carma_obj_handler = (caObjZ *) (handler->carma_object);
			delete carma_obj_handler;
		} else if (handler->type == Y_SCOMPLEX) {
			caObjC *carma_obj_handler = (caObjC *) (handler->carma_object);
			delete carma_obj_handler;
		} else
			throw "Type unknown";
		handler->type = Y_VOID;
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	}
}

void Y_yoga_obj(int argc)
/** @brief yObj_struct constructor.
 *  @param[in] argc : command line argument(s)
 *  several cases are handled :
 *   - a pointer to a yoga_obj is passed as an argument : copy the pointee in a new yoga_obj
 *   - a yoga_obj is passed as an argument : create a copy
 *   - a string (type of object) and an array of dimensions (Yorick convention) is passed
 *   - an existing array is passed : create a new object an fill it with array content
 *  supported types : float, double, scomplex, complex, uint, int, long
 *  the created object is pushed on the stack
 */
{
	long ntot;
	long dims[Y_DIMSIZE];
	void *data_input = 0L;

	try {
		int oType;
		int yType = yarg_typeid(argc - 1);
		if (yType == Y_POINTER) {
			yObj_struct *handle_obj = (yObj_struct *) ygets_p(argc - 1);
			yObj_struct *handle = (yObj_struct *) ypush_obj(&yObj,
					sizeof(yObj_struct));
			handle->type = handle_obj->type;
			handle->device = handle_obj->device;
			handle->isRef = 1;
			handle->carma_object = handle_obj->carma_object;
		} else if (yType == Y_OPAQUE) { // Copy constructor
			yObj_struct *handle_obj = (yObj_struct *) yget_obj(argc - 1,
					&yObj);
			yObj_struct *handle = (yObj_struct *) ypush_obj(&yObj,
					sizeof(yObj_struct));
			handle->type = handle_obj->type;
			handle->device = handle_obj->device;
			handle->isRef = 0;
			carma_context *context_handle = _getCurrentContext();
			context_handle->set_activeDevice(handle->device);
			if (handle->type == Y_FLOAT) {
				caObjS *carma_obj_handler = (caObjS *) (handle_obj->carma_object);
				handle->carma_object = new caObjS(context_handle,
						carma_obj_handler);
			} else if (handle->type == Y_DOUBLE) {
				caObjD *carma_obj_handler = (caObjD *) (handle_obj->carma_object);
				handle->carma_object = new caObjD(context_handle,
						carma_obj_handler);
			} else if (handle->type == Y_INT) {
				caObjI *carma_obj_handler = (caObjI *) (handle_obj->carma_object);
				handle->carma_object = new caObjI(context_handle,
						carma_obj_handler);
			} else if (handle->type == Y_COMPLEX) {
				caObjZ *carma_obj_handler = (caObjZ *) (handle_obj->carma_object);
				handle->carma_object = new caObjZ(context_handle,
						carma_obj_handler);
			} else if (handle->type == Y_SCOMPLEX) {
				caObjC *carma_obj_handler = (caObjC *) (handle_obj->carma_object);
				handle->carma_object = new caObjC(context_handle,
						carma_obj_handler);
			}
		} else { // Standard constructor
			carma_context *context_handle = _getCurrentContext();
			int activeDevice = context_handle->get_activeDevice();
			int odevice = activeDevice;
			if (yType == Y_STRING) { // based on the input type and dimensions array
				data_input = ygeta_l(argc - 2, &ntot, dims);
				char *type_data = ygets_q(argc - 1);
				if (strcmp(type_data, "float") == 0) {
					oType = Y_FLOAT;
				} else if (strcmp(type_data, "double") == 0) {
					oType = Y_DOUBLE;
				} else if (strcmp(type_data, "int") == 0) {
					oType = Y_INT;
				} else if (strcmp(type_data, "complex") == 0) {
					oType = Y_COMPLEX;
				} else if (strcmp(type_data, "scomplex") == 0) {
					oType = Y_SCOMPLEX;
				} else if (strcmp(type_data, "uint") == 0) {
					oType = Y_SHORT;
				}
				if (argc > 2)
					odevice = ygets_i(argc-3);
			} else { // based on the input array
				data_input = ygeta_any(argc - 1, &ntot, dims, &oType);
				if (argc > 1)
					odevice = ygets_i(argc-2);
			}

			activeDevice = context_handle->set_activeDevice(odevice);

			yObj_struct *handle = (yObj_struct *) ypush_obj(&yObj,
					sizeof(yObj_struct));
			handle->device = odevice;
			handle->type = oType;
			handle->isRef = 0;

			if (oType == Y_FLOAT) {
				if (yType == Y_STRING) { // based on the input type and dimensions array
					handle->carma_object = new caObjS(context_handle,
							(long*) data_input);
				} else { // based on the input array
					handle->carma_object = new caObjS(context_handle, dims,
							(float*) data_input);
				}
			} else if (oType == Y_DOUBLE) {
				if (yType == Y_STRING) {
					handle->carma_object = new caObjD(context_handle,
							(long*) data_input);
				} else {
					handle->carma_object = new caObjD(context_handle, dims,
							(double*) data_input);
				}
			} else if (oType == Y_INT) {
				if (yType == Y_STRING) {
					handle->carma_object = new caObjI(context_handle,
							(long*) data_input);
				} else {
					handle->carma_object = new caObjI(context_handle, dims,
							(int*) data_input);
				}
			} else if (oType == Y_SHORT) {
				if (yType == Y_STRING) {
					handle->carma_object = new caObjUI(context_handle,
							(long*) data_input);
				} else {
					handle->carma_object = new caObjUI(context_handle, dims,
							(uint*) data_input);
				}
			} else if (oType == Y_COMPLEX) {
				if (yType == Y_STRING) {
					handle->carma_object = new caObjZ(context_handle,
							(long*) data_input);
				} else {
					handle->carma_object = new caObjZ(context_handle, dims,
							(cuDoubleComplex*) data_input);
				}
			} else if (oType == Y_SCOMPLEX) {
				if (yType == Y_STRING) {
					handle->carma_object = new caObjC(context_handle,
							(long*) data_input);
				} else {
					stringstream buf;
					buf << "Type not supported in yorick in " << __FILE__ << "@"
							<< __LINE__ << endl;
					throw buf.str();
					//handle->carma_object = new caObjC( dims, (cuFloatComplex*)data_input);
				}
			} else {
				stringstream buf;
				buf << "Type not found in " << __FILE__ << "@" << __LINE__
						<< endl;
				throw buf.str();
			}
		}
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	} catch (...) {
		stringstream buf;
		buf << "unknown error with yObj construction in " << __FILE__ << "@"
				<< __LINE__ << endl;
		y_error(buf.str().c_str());
	}
}

/*
 *         _   _ _ _ _   _
 *   _   _| |_(_) (_) |_(_) ___  ___
 *  | | | | __| | | | __| |/ _ \/ __|
 *  | |_| | |_| | | | |_| |  __/\__ \
 *   \__,_|\__|_|_|_|\__|_|\___||___/
 */

void Y_yoga_getp(int argc)
/** @brief simple routine to retrieve the address to a yoga_obj
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  the corresponding address is pushed on the stack
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	long a = 0;
	ypointer_t *ptr_handle = ypush_p(&a);
	*ptr_handle = handle;
}

void Y_yoga_getpl(int argc)
/** @brief simple routine to retrieve the address to a yoga_obj
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  the corresponding address is pushed on the stack
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	//long a = 0;
	ypush_long((long)handle);
}



void Y_yoga_setv(int argc)
/** @brief simple routine to fill a vector yoga_obj with data
 *  @param[in] argc : command line arguments
 *    - first  : a yoga_obj
 *    - second : the data (see yObj_struct constructor for supported types)
 */
{
	long ntot;
	long dims[Y_DIMSIZE];

	void *data_input = 0L;
	int yType = yarg_typeid(argc - 1);

	data_input = ygeta_any(0, &ntot, dims, &yType);

	yObj_struct *handle = (yObj_struct *) ypush_obj(&yObj,
			sizeof(yObj_struct));
	carma_context *context_handle = _getCurrentContext();

	context_handle->set_activeDevice(handle->device);
	handle->type = yType;

	long dims_data[2];
	dims_data[0] = 1;
	if (dims[0] == 1)
		dims_data[1] = dims[1];
	if (dims[0] > 1)
		dims_data[1] = ntot;

	if (yType == Y_FLOAT) {
		handle->carma_object = new caObjS(context_handle, (long*) dims_data);
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		carma_obj_handler->host2deviceVect((float *) data_input, 1, 1);
	} else if (yType == Y_DOUBLE) {
		handle->carma_object = new caObjD(context_handle, (long*) dims_data);
		caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
		carma_obj_handler->host2deviceVect((double *) data_input, 1, 1);
	} else {
		y_error("Type not supported\n");
	}
}

void Y_yoga_setm(int argc)
/** @brief simple routine to fill a matrix yoga_obj with data
 *  @param[in] argc : command line arguments
 *    - first  : a yoga_obj
 *    - second : the data (see yObj_struct constructor for supported types)
 */
{
	long ntot;
	long dims[Y_DIMSIZE];

	void *data_input = 0L;
	int yType = yarg_typeid(argc - 1);
	data_input = ygeta_any(0, &ntot, dims, &yType);

	yObj_struct *handle = (yObj_struct *) ypush_obj(&yObj,
			sizeof(yObj_struct));

	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(handle->device);
	handle->type = yType;
	if (yType == Y_FLOAT) {
		handle->carma_object = new caObjS(context_handle, (long*) dims);
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		carma_obj_handler->host2deviceMat((float *) data_input, dims[1],
				dims[1]);
	} else if (yType == Y_DOUBLE) {
		handle->carma_object = new caObjD(context_handle, (long*) dims);
		caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
		carma_obj_handler->host2deviceMat((double *) data_input, dims[1],
				dims[1]);
	} else {
		y_error("Type not supported\n");
	}
}

void Y_yoga_host2device(int argc)
/** @brief simple routine to fill a yoga_obj with data
 *  @param[in] argc : command line arguments
 *    - first  : a yoga_obj
 *    - second : the data (see yObj_struct constructor for supported types)
 */
{
	if (yarg_subroutine()) {
		try {
			yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1,
					&yObj);
			carma_context *context_handle = _getCurrentContext();
			context_handle->set_activeDeviceForCpy(handle->device);
			long ntot;
			long dims;
			if (handle->type == Y_FLOAT) {
				caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
				float *data = ygeta_f(argc - 2, &ntot, &dims);
				carma_obj_handler->host2device(data);
			} else if (handle->type == Y_DOUBLE) {
				caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
				double *data = ygeta_d(argc - 2, &ntot, &dims);
				carma_obj_handler->host2device(data);
			} else if (handle->type == Y_INT) {
				caObjI *carma_obj_handler = (caObjI *) (handle->carma_object);
				int *data = ygeta_i(argc - 2, &ntot, &dims);
				carma_obj_handler->host2device(data);
			} else if (handle->type == Y_SHORT) {
				caObjUI *carma_obj_handler = (caObjUI *) (handle->carma_object);
				unsigned int *data = (unsigned int *) ygeta_i(argc - 2, &ntot,
						&dims);
				carma_obj_handler->host2device(data);
			} else if (handle->type == Y_SCOMPLEX) {
				caObjC *carma_obj_handler = (caObjC *) (handle->carma_object);
				float *data = ygeta_f(argc - 2, &ntot, &dims);
				carma_obj_handler->host2device((cuFloatComplex*) data);
			} else if (handle->type == Y_COMPLEX) {
				caObjZ *carma_obj_handler = (caObjZ *) (handle->carma_object);
				double *data = ygeta_d(argc - 2, &ntot, &dims);
				carma_obj_handler->host2device((cuDoubleComplex*) data);
			} else
				throw "Type not found";
		} catch (string &msg) {
			y_error(msg.c_str());
		} catch (char const * msg) {
			y_error(msg);
		} catch (...) {
			y_error("unknown error in Y_yoga_host2device\n");
		}

	} else {
		y_error(
				"yoga_host2device only as a subroutine \n use yoga_obj() if you want to create a gpu object from a yorick array");
	}
}

void Y_yoga_device2host(int argc)
/** @brief simple routine to retrieve data from a yoga_obj
 *  @param[in] argc : command line arguments
 *    - first  : a yoga_obj
 *    - second : the data (see yObj_struct constructor for supported types)
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDeviceForCpy(handle->device);
	int opt;
	if (argc > 1)
		opt = ygets_i(argc-2);
	else
		opt = 0;
	if ((opt != 1) && (opt != 0))
		y_error("expecting optioal flag 1/0");

	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		float *data = ypush_f(carma_obj_handler->getDims());
		if (opt == 0)
			carma_obj_handler->device2host(data);
		else
			carma_obj_handler->device2hostOpt(data);
	} else if (handle->type == Y_DOUBLE) {
		caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
		double *data = ypush_d((long *) carma_obj_handler->getDims());
		if (opt == 0)
			carma_obj_handler->device2host(data);
		else
			carma_obj_handler->device2hostOpt(data);
	}
}

/*
 *  _     _
 * | |__ | | __ _ ___
 * | '_ \| |/ _` / __|
 * | |_) | | (_| \__ \
 * |_.__/|_|\__,_|___/
 */

void Y_yoga_imin(int argc)
/** @brief wrapper routine for yoga_blas imin method
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  push the result on the stack
 *  only floating point types supported (single or double precision)
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(handle->device);
	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		ypush_int(carma_obj_handler->imin(1)); // here 1 is the increment
	} else if (handle->type == Y_DOUBLE) {
		caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
		ypush_int(carma_obj_handler->imin(1)); // here 1 is the increment
	}
}

void Y_yoga_imax(int argc)
/** @brief wrapper routine for yoga_blas imax method
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  push the result on the stack
 *  only floating point types supported (single or double precision)
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(handle->device);

	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		ypush_int(carma_obj_handler->imax(1)); // here 1 is the increment
	}
	if (handle->type == Y_DOUBLE) {
		caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
		ypush_int(carma_obj_handler->imax(1)); // here 1 is the increment
	}
}

void Y_yoga_asum(int argc)
/** @brief wrapper routine for yoga_blas asum method
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  push the result on the stack
 *  only floating point types supported (single or double precision)
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDeviceForCpy(handle->device);
	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		ypush_double((double) carma_obj_handler->asum(1)); // here 1 is the increment
	}
	if (handle->type == Y_DOUBLE) {
		caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
		ypush_double(carma_obj_handler->asum(1)); // here 1 is the increment
	}
}

void Y_yoga_sum(int argc)
/** @brief wrapper routine for yoga sum method
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  push the result on the stack
 *  only floating point types supported (single or double precision)
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDeviceForCpy(handle->device);
	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		ypush_double((double) carma_obj_handler->sum()); 
	}
	if (handle->type == Y_DOUBLE) {
		caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
		ypush_double(carma_obj_handler->sum()); 
	}
}

void Y_yoga_nrm2(int argc)
/** @brief wrapper routine for yoga_blas nrm2 method
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  push the result on the stack
 *  only floating point types supported (single or double precision)
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(handle->device);
	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		ypush_double((double) carma_obj_handler->nrm2(1)); // here 1 is the increment
	}
	if (handle->type == Y_DOUBLE) {
		caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
		ypush_double(carma_obj_handler->nrm2(1)); // here 1 is the increment
	}
}

void Y_yoga_scale(int argc)
/** @brief wrapper routine for yoga_blas scale method
 *  @param[in] argc : command line arguments
 *    - first  : a yoga_obj
 *    - second : the scaling factor
 *  only floating point types supported (single or double precision)
 *  inplace only
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDeviceForCpy(handle->device);
		if (handle->type == Y_FLOAT) {
			float alpha = ygets_f(argc-2);
			caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
			carma_obj_handler->scale(alpha, 1);
		}
		if (handle->type == Y_DOUBLE) {
			double alpha = ygets_d(argc - 2);
			caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
			carma_obj_handler->scale(alpha, 1);
		}
	} else {
		y_error("yoga_scale only inplace (subroutine) \n");
	}
}

void Y_yoga_swap(int argc)
/** @brief wrapper routine for yoga_blas swap method
 *  @param[in] argc : command line arguments
 *    - first  : a yoga_obj
 *    - second : another yoga_obj
 *  only floating point types supported (single or double precision)
 */
{
	yObj_struct *handle_dest = (yObj_struct *) yget_obj(argc - 1, &yObj);
	yObj_struct *handle_src = (yObj_struct *) yget_obj(argc - 2, &yObj);
	if (handle_dest->device != handle_src->device)
		y_error("swap only on the same device");
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDeviceForCpy(handle_dest->device);
	if (handle_dest->type == Y_FLOAT) {
		caObjS *carma_obj_handler_dest = (caObjS *) (handle_dest->carma_object);
		caObjS *carma_obj_handler_src = (caObjS *) (handle_src->carma_object);
		carma_obj_handler_dest->swap(carma_obj_handler_src, 1, 1);
	}
	if (handle_dest->type == Y_DOUBLE) {
		caObjD *carma_obj_handler_dest = (caObjD *) (handle_dest->carma_object);
		caObjD *carma_obj_handler_src = (caObjD *) (handle_src->carma_object);
		carma_obj_handler_dest->swap(carma_obj_handler_src, 1, 1);
	}
}

void Y_yoga_axpy(int argc)
/** @brief wrapper routine for yoga_blas axpy method
 *  @param[in] argc : command line arguments
 *  can work as a (1) subroutine (return discarded) or (2) as a function
 *    - first  : (1) the destnation yoga_obj / (2) the source yoga_obj
 *    - second : (1) / (2) scaling factor
 *    - third  : (1) the source yoga_obj
 *  only floating point types supported (single or double precision)
 */
{
	if (yarg_subroutine()) {
		// called as a subroutine : destination already exists
		yObj_struct *handle_dest = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_src = (yObj_struct *) yget_obj(argc - 3,
				&yObj);
		if (handle_dest->device != handle_src->device)
			y_error("axpy only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDeviceForCpy(handle_dest->device);
		if (handle_src->type == Y_FLOAT) {
			float alpha = ygets_f(argc-2);
			caObjS *carma_obj_handler_dest = (caObjS *) (handle_dest->carma_object);
			caObjS *carma_obj_handler_src = (caObjS *) (handle_src->carma_object);
			carma_obj_handler_dest->axpy(alpha, carma_obj_handler_src, 1, 1);
		} else if (handle_src->type == Y_DOUBLE) {
			double alpha = ygets_d(argc - 2);
			caObjD *carma_obj_handler_dest = (caObjD *) (handle_dest->carma_object);
			caObjD *carma_obj_handler_src = (caObjD *) (handle_src->carma_object);
			carma_obj_handler_dest->axpy(alpha, carma_obj_handler_src, 1, 1);
		}
	} else {
		// called as a function : we need to create a destination object
		yObj_struct *handle_src = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDeviceForCpy(handle_src->device);
		yObj_struct *handle_dest = (yObj_struct *) ypush_obj(&yObj,
				sizeof(yObj_struct));
		handle_dest->device = handle_src->device;
		if (handle_src->type == Y_FLOAT) {
			float alpha = ygets_f(argc-2);
			handle_dest->type = handle_src->type;
			caObjS *carma_obj_handler_src = (caObjS *) (handle_src->carma_object);
			handle_dest->carma_object = new caObjS(context_handle,
					carma_obj_handler_src->getDims());
			caObjS *carma_obj_handler_dest = (caObjS *) (handle_dest->carma_object);
			carma_obj_handler_dest->axpy(alpha, carma_obj_handler_src, 1, 1);
		} else if (handle_src->type == Y_DOUBLE) {
			double alpha = ygets_d(argc - 2);
			handle_dest->type = handle_src->type;
			caObjD *carma_obj_handler_src = (caObjD *) (handle_src->carma_object);
			handle_dest->carma_object = new caObjD(context_handle,
					carma_obj_handler_src->getDims());
			caObjD *carma_obj_handler_dest = (caObjD *) (handle_dest->carma_object);
			carma_obj_handler_dest->axpy(alpha, carma_obj_handler_src, 1, 1);
		}
	}
}

void Y_yoga_dot(int argc)
/** @brief wrapper routine for yoga_blas dot method
 *  @param[in] argc : command line arguments
 *    - first  : a yoga_obj
 *    - second : another yoga_obj
 *  push the result on the stack
 *  only floating point types supported (single or double precision)
 */
{
	yObj_struct *handle1 = (yObj_struct *) yget_obj(argc - 1, &yObj);
	yObj_struct *handle2 = (yObj_struct *) yget_obj(argc - 2, &yObj);
	if (handle1->device != handle2->device)
		y_error("dot only on the same device");
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDeviceForCpy(handle2->device);
	if (handle1->type == Y_FLOAT) {
		caObjS *carma_obj_handler1 = (caObjS *) (handle1->carma_object);
		caObjS *carma_obj_handler2 = (caObjS *) (handle2->carma_object);
		ypush_double(carma_obj_handler1->dot(carma_obj_handler2, 1, 1));
		// here 1 is the increment
	}
	if (handle1->type == Y_DOUBLE) {
		caObjD *carma_obj_handler1 = (caObjD *) (handle1->carma_object);
		caObjD *carma_obj_handler2 = (caObjD *) (handle2->carma_object);
		ypush_double(carma_obj_handler1->dot(carma_obj_handler2, 1, 1));
		// here 1 is the increment
	}
}

void Y_yoga_mv(int argc)
/** @brief wrapper routine for yoga_blas mv method
 *  @param[in] argc : command line arguments
 *  can work as a (1) subroutine (return discarded) or (2) as a function
 *    - first  : (1) the destnation vector yoga_obj / (2) the matrix yoga_obj
 *    - second : (1) the matrix yoga_obj / (2) the source vector yoga_obj
 *    - third  : (1) the source vector yoga_obj / (2) optional scaling factor for dest
 *    - fourth : (1) optional scaling factor for dest / (2) optional scaling factor for src
 *    - fifth  : (1) optional scaling factor for src
 *  in case (2) the destination is pushed on the stack as a yoga_obj
 *  only floating point types supported (single or double precision)
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle_vecty = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_mat = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		yObj_struct *handle_vectx = (yObj_struct *) yget_obj(argc - 3,
				&yObj);
		if ((handle_vecty->device != handle_mat->device)
				|| (handle_vecty->device != handle_vectx->device)
				|| (handle_mat->device != handle_vectx->device))
			y_error("mv only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_mat->device);

		if (handle_mat->type == Y_FLOAT) {
			caObjS *carma_obj_handler_mat = (caObjS *) (handle_mat->carma_object);
			caObjS *carma_obj_handler_vectx =
					(caObjS *) (handle_vectx->carma_object);
			caObjS *carma_obj_handler_vecty =
					(caObjS *) (handle_vecty->carma_object);
			float alpha;
			if (argc > 3) {
				alpha = ygets_f(argc-4);
			} else
				alpha = 1.0f;
			float beta;
			if (argc > 4) {
				beta = ygets_f(argc-5);
			} else
				beta = 0.0f;
			carma_obj_handler_vecty->gemv('n', alpha, carma_obj_handler_mat,
					carma_obj_handler_mat->getDims(1), carma_obj_handler_vectx, 1,
					beta, 1);
			// here 1 is the increment
		}
		if (handle_mat->type == Y_DOUBLE) {
			caObjD *carma_obj_handler_mat = (caObjD *) (handle_mat->carma_object);
			caObjD *carma_obj_handler_vectx =
					(caObjD *) (handle_vectx->carma_object);
			caObjD *carma_obj_handler_vecty =
					(caObjD *) (handle_vecty->carma_object);
			double alpha;
			if (argc > 3) {
				alpha = ygets_d(argc - 4);
			} else
				alpha = 1.0;
			double beta;
			if (argc > 4) {
				beta = ygets_d(argc - 5);
			} else
				beta = 0.0;
			carma_obj_handler_vecty->gemv('n', alpha, carma_obj_handler_mat,
					carma_obj_handler_mat->getDims(1), carma_obj_handler_vectx, 1,
					beta, 1);
			// here 1 is the increment
		}
	} else {
		// called as a function : need to allocate space
		yObj_struct *handle_mat = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_vectx = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		if (handle_vectx->device != handle_mat->device)
			y_error("mv only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_mat->device);
		yObj_struct *handle_vecty = (yObj_struct *) ypush_obj(&yObj,
				sizeof(yObj_struct));
		handle_vecty->device = handle_mat->device;
		if (handle_mat->type == Y_FLOAT) {
			float alpha;
			if (argc > 2) {
				alpha = ygets_f(argc-3);
			} else
				alpha = 1.0f;
			float beta;
			if (argc > 3) {
				beta = ygets_f(argc-4);
			} else
				beta = 0.0f;
			handle_vecty->type = handle_vectx->type;
			caObjS *carma_obj_handler_mat = (caObjS *) (handle_mat->carma_object);
			caObjS *carma_obj_handler_vectx =
					(caObjS *) (handle_vectx->carma_object);
			long dims_data_y[2];
			dims_data_y[0] = 1;
			dims_data_y[1] = carma_obj_handler_mat->getDims(1);
			handle_vecty->carma_object = new caObjS(context_handle, dims_data_y);
			caObjS *carma_obj_handler_vecty =
					(caObjS *) (handle_vecty->carma_object);
			carma_obj_handler_vecty->gemv('n', 1.0f, carma_obj_handler_mat,
					carma_obj_handler_mat->getDims(1), carma_obj_handler_vectx, 1,
					0.0f, 1);
			// here 1 is the increment
		} else if (handle_mat->type == Y_DOUBLE) {
			handle_vecty->type = handle_vectx->type;
			caObjD *carma_obj_handler_mat = (caObjD *) (handle_mat->carma_object);
			caObjD *carma_obj_handler_vectx =
					(caObjD *) (handle_vectx->carma_object);
			double alpha;
			if (argc > 2) {
				alpha = ygets_d(argc - 3);
			} else
				alpha = 1.0;
			double beta;
			if (argc > 3) {
				beta = ygets_d(argc - 4);
			} else
				beta = 0.0;
			long dims_data_y[2];
			dims_data_y[0] = 1;
			dims_data_y[1] = carma_obj_handler_mat->getDims(1);
			handle_vecty->carma_object = new caObjD(context_handle, dims_data_y);
			caObjD *carma_obj_handler_vecty =
					(caObjD *) (handle_vecty->carma_object);
			carma_obj_handler_vecty->gemv('n', 1.0, carma_obj_handler_mat,
					carma_obj_handler_mat->getDims(1), carma_obj_handler_vectx, 1,
					0.0, 1);
			// here 1 is the increment
		}
	}
}

void Y_yoga_rank1(int argc)
/** @brief wrapper routine for yoga_blas rank1 method
 *  @param[in] argc : command line arguments
 *  can work as a (1) subroutine (return discarded) or (2) as a function
 *    - first  : (1) the destination matrix yoga_obj / (2) the x vector yoga_obj
 *    - second : (1) the x vector yoga_obj / (2) the y vector yoga_obj
 *    - third  : (1) the y vector yoga_obj
 *  in case (2) the destination is pushed on the stack as a yoga_obj
 *  only floating point types supported (single or double precision)
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle_mat = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_vectx = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		yObj_struct *handle_vecty = (yObj_struct *) yget_obj(argc - 3,
				&yObj);
		if ((handle_vecty->device != handle_mat->device)
				|| (handle_vecty->device != handle_vectx->device)
				|| (handle_mat->device != handle_vectx->device))
			y_error("rank1 only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_mat->device);
		if (handle_mat->type == Y_FLOAT) {
			caObjS *carma_obj_handler_mat = (caObjS *) (handle_mat->carma_object);
			caObjS *carma_obj_handler_vectx =
					(caObjS *) (handle_vectx->carma_object);
			caObjS *carma_obj_handler_vecty =
					(caObjS *) (handle_vecty->carma_object);
			carma_obj_handler_mat->ger(1.0f, carma_obj_handler_vectx, 1,
					carma_obj_handler_vecty, 1,
					carma_obj_handler_mat->getDims(1));
			// here 1 is the increment
		}
		if (handle_mat->type == Y_DOUBLE) {
			caObjD *carma_obj_handler_mat = (caObjD *) (handle_mat->carma_object);
			caObjD *carma_obj_handler_vectx =
					(caObjD *) (handle_vectx->carma_object);
			caObjD *carma_obj_handler_vecty =
					(caObjD *) (handle_vecty->carma_object);
			carma_obj_handler_mat->ger(1.0, carma_obj_handler_vectx, 1,
					carma_obj_handler_vecty, 1,
					carma_obj_handler_mat->getDims(1));
			// here 1 is the increment
		}
	} else {
		// called as a function : need to allocate space
		yObj_struct *handle_vectx = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_vecty = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		if (handle_vectx->device != handle_vecty->device)
			y_error("rank1 only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_vectx->device);
		yObj_struct *handle_matA = (yObj_struct *) ypush_obj(&yObj,
				sizeof(yObj_struct));
		handle_matA->device = handle_vecty->device;
		if (handle_vectx->type == Y_FLOAT) {
			handle_matA->type = handle_vectx->type;
			caObjS *carma_obj_handler_vecty =
					(caObjS *) (handle_vecty->carma_object);
			caObjS *carma_obj_handler_vectx =
					(caObjS *) (handle_vectx->carma_object);
			long dims_data_mat[3];
			dims_data_mat[0] = 1;
			dims_data_mat[1] = carma_obj_handler_vecty->getDims(1);
			dims_data_mat[2] = carma_obj_handler_vectx->getDims(1);
			handle_matA->carma_object = new caObjS(context_handle, dims_data_mat);
			caObjS *carma_obj_handler_matA = (caObjS *) (handle_matA->carma_object);
			carma_obj_handler_matA->ger(1.0f, carma_obj_handler_vectx, 1,
					carma_obj_handler_vecty, 1,
					carma_obj_handler_matA->getDims(1));
			// here 1 is the increment
		} else if (handle_vectx->type == Y_DOUBLE) {
			handle_matA->type = handle_vectx->type;
			caObjD *carma_obj_handler_vecty =
					(caObjD *) (handle_vecty->carma_object);
			caObjD *carma_obj_handler_vectx =
					(caObjD *) (handle_vectx->carma_object);
			long dims_data_mat[3];
			dims_data_mat[0] = 1;
			dims_data_mat[1] = carma_obj_handler_vecty->getDims(1);
			dims_data_mat[2] = carma_obj_handler_vectx->getDims(1);
			handle_matA->carma_object = new caObjD(context_handle, dims_data_mat);
			caObjD *carma_obj_handler_matA = (caObjD *) (handle_matA->carma_object);
			carma_obj_handler_matA->ger(1.0, carma_obj_handler_vectx, 1,
					carma_obj_handler_vecty, 1,
					carma_obj_handler_matA->getDims(1));
			// here 1 is the increment
		}
	}
}

void Y_yoga_mm(int argc)
/** @brief wrapper routine for yoga_blas mm method
 *  @param[in] argc : command line arguments
 *  can work as a (1) subroutine (return discarded) or (2) as a function
 *    - first   : (1) the C matrix yoga_obj / (2) the A matrix yoga_obj
 *    - second  : (1) the A matrix yoga_obj / (2) the B matrix yoga_obj
 *    - third   : (1) the B matrix yoga_obj
 TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *    - fourth  : (1) the alpha coeff
 *    - fifth   : (1) the beta coeff
 *    - sixth   : (1) the opA
 *    - seventh : (1) the opB
 *  in case (2) the destination is pushed on the stack as a yoga_obj
 *  only floating point types supported (single or double precision)
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle_matC = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_matA = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		yObj_struct *handle_matB = (yObj_struct *) yget_obj(argc - 3,
				&yObj);
		if ((handle_matA->device != handle_matB->device)
				|| (handle_matA->device != handle_matC->device)
				|| (handle_matB->device != handle_matC->device))
			y_error("mm only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_matA->device);

		char opA = 'n';
		if (argc > 3)
			opA = ygets_c(argc-4);
		char opB = 'n';
		if (argc > 4)
			opB = ygets_c(argc-5);

		if (handle_matC->type == Y_FLOAT) {
			float alpha = 1.0f;
			if (argc > 5)
				alpha = ygets_f(argc-6);
			float beta = 0.0f;
			if (argc > 6)
				beta = ygets_f(argc-7);

			caObjS *carma_obj_handler_matA = (caObjS *) (handle_matA->carma_object);
			caObjS *carma_obj_handler_matB = (caObjS *) (handle_matB->carma_object);
			caObjS *carma_obj_handler_matC = (caObjS *) (handle_matC->carma_object);

			carma_obj_handler_matC->gemm(opA, opB, alpha, carma_obj_handler_matA,
					carma_obj_handler_matA->getDims(1), carma_obj_handler_matB,
					carma_obj_handler_matB->getDims(1), beta,
					carma_obj_handler_matC->getDims(1));
		}
		if (handle_matC->type == Y_DOUBLE) {
			double alpha = 1.0;
			if (argc > 5)
				alpha = ygets_d(argc - 6);
			double beta = 0.0;
			if (argc > 6)
				beta = ygets_d(argc - 7);

			caObjD *carma_obj_handler_matA = (caObjD *) (handle_matA->carma_object);
			caObjD *carma_obj_handler_matB = (caObjD *) (handle_matB->carma_object);
			caObjD *carma_obj_handler_matC = (caObjD *) (handle_matC->carma_object);

			carma_obj_handler_matC->gemm(opA, opB, alpha, carma_obj_handler_matA,
					carma_obj_handler_matA->getDims(1), carma_obj_handler_matB,
					carma_obj_handler_matB->getDims(1), beta,
					carma_obj_handler_matC->getDims(1));
		}
	} else {
		// called as a function : need to allocate space
		yObj_struct *handle_matA = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_matB = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		if (handle_matA->device != handle_matB->device)
			y_error("mm only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_matA->device);

		char opA = 'n';
		if (argc > 2)
			opA = ygets_c(argc-3);
		char opB = 'n';
		if (argc > 3)
			opB = ygets_c(argc-4);

		if (handle_matA->type == Y_FLOAT) {
			float alpha = 1.0f;
			if (argc > 4)
				alpha = ygets_f(argc-5);
			float beta = 0.0f;
			if (argc > 5)
				beta = ygets_f(argc-6);

			yObj_struct *handle_matC = (yObj_struct *) ypush_obj(&yObj,
					sizeof(yObj_struct));
			handle_matC->device = handle_matA->device;

			handle_matC->type = handle_matA->type;
			caObjS *carma_obj_handler_matA = (caObjS *) (handle_matA->carma_object);
			caObjS *carma_obj_handler_matB = (caObjS *) (handle_matB->carma_object);

			long dims_data_mat[3];
			dims_data_mat[0] = 2;
			if (opA == 'n')
				dims_data_mat[1] = carma_obj_handler_matA->getDims(1);
			else
				dims_data_mat[1] = carma_obj_handler_matA->getDims(2);
			if (opB == 'n')
				dims_data_mat[2] = carma_obj_handler_matB->getDims(2);
			else
				dims_data_mat[2] = carma_obj_handler_matB->getDims(1);

			handle_matC->carma_object = new caObjS(context_handle, dims_data_mat);
			caObjS *carma_obj_handler_matC = (caObjS *) (handle_matC->carma_object);

			carma_obj_handler_matC->gemm(opA, opB, alpha, carma_obj_handler_matA,
					carma_obj_handler_matA->getDims(1), carma_obj_handler_matB,
					carma_obj_handler_matB->getDims(1), beta,
					carma_obj_handler_matC->getDims(1));
		} else if (handle_matA->type == Y_DOUBLE) {
			double alpha = 1.0;
			if (argc > 4)
				alpha = ygets_d(argc - 5);
			double beta = 0.0;
			if (argc > 5)
				beta = ygets_d(argc - 6);

			yObj_struct *handle_matC = (yObj_struct *) ypush_obj(&yObj,
					sizeof(yObj_struct));
			handle_matC->device = handle_matA->device;

			handle_matC->type = handle_matA->type;
			caObjD *carma_obj_handler_matA = (caObjD *) (handle_matA->carma_object);
			caObjD *carma_obj_handler_matB = (caObjD *) (handle_matB->carma_object);

			long dims_data_mat[3];
			dims_data_mat[0] = 2;
			if (opA == 'n')
				dims_data_mat[1] = carma_obj_handler_matA->getDims(1);
			else
				dims_data_mat[1] = carma_obj_handler_matA->getDims(2);
			if (opB == 'n')
				dims_data_mat[2] = carma_obj_handler_matB->getDims(2);
			else
				dims_data_mat[2] = carma_obj_handler_matB->getDims(1);

			handle_matC->carma_object = new caObjD(context_handle, dims_data_mat);
			caObjD *carma_obj_handler_matC = (caObjD *) (handle_matC->carma_object);
			carma_obj_handler_matC->gemm(opA, opB, alpha, carma_obj_handler_matA,
					carma_obj_handler_matA->getDims(1), carma_obj_handler_matB,
					carma_obj_handler_matB->getDims(1), beta,
					carma_obj_handler_matC->getDims(1));
		}
	}
}

/*
 *   _
 *  | |_ _ __ __ _ _ __  ___ _ __   ___  ___  ___
 *  | __| '__/ _` | '_ \/ __| '_ \ / _ \/ __|/ _ \
 *  | |_| | | (_| | | | \__ \ |_) | (_) \__ \  __/
 *   \__|_|  \__,_|_| |_|___/ .__/ \___/|___/\___|
 *                          |_|
 */

void Y_yoga_transpose(int argc)
/** @brief wrapper routine for yoga transpose method
 *  @param[in] argc : command line arguments
 *  can work as a (1) subroutine (return discarded) or (2) as a function
 *    - first  : (1) the destination matrix yoga_obj / (2) the A matrix yoga_obj
 *    - second : (1) the source matrix yoga_obj
 *  in case (2) the destination is pushed on the stack as a yoga_obj
 *  only floating point types supported (single or double precision)
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle_dest = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_src = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		if (handle_dest->device != handle_src->device)
			y_error("transpose only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDeviceForCpy(handle_dest->device);
		if (handle_src->type == Y_FLOAT) {
			caObjS *carma_obj_handler_dest = (caObjS *) (handle_dest->carma_object);
			caObjS *carma_obj_handler_src = (caObjS *) (handle_src->carma_object);
			carma_obj_handler_dest->transpose(carma_obj_handler_src);
		} else if (handle_src->type == Y_DOUBLE) {
			caObjD *carma_obj_handler_dest = (caObjD *) (handle_dest->carma_object);
			caObjD *carma_obj_handler_src = (caObjD *) (handle_src->carma_object);
			carma_obj_handler_dest->transpose(carma_obj_handler_src);
		}
	} else {
		// called as a function : we need to create an object
		yObj_struct *handle_src = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_dest = (yObj_struct *) ypush_obj(&yObj,
				sizeof(yObj_struct));
		handle_dest->device = handle_src->device;
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDeviceForCpy(handle_dest->device);
		if (handle_src->type == Y_FLOAT) {
			handle_dest->type = handle_src->type;
			caObjS *carma_obj_handler_src = (caObjS *) (handle_src->carma_object);
			long dims_data_dest[3];
			dims_data_dest[0] = carma_obj_handler_src->getDims(0);
			dims_data_dest[1] = carma_obj_handler_src->getDims(2);
			dims_data_dest[2] = carma_obj_handler_src->getDims(1);
			handle_dest->carma_object = new caObjS(context_handle,
					dims_data_dest);
			caObjS *carma_obj_handler_dest = (caObjS *) (handle_dest->carma_object);
			carma_obj_handler_dest->transpose(carma_obj_handler_src);
		} else if (handle_src->type == Y_DOUBLE) {
			handle_dest->type = handle_src->type;
			caObjD *carma_obj_handler_src = (caObjD *) (handle_src->carma_object);
			long dims_data_dest[3];
			dims_data_dest[0] = carma_obj_handler_src->getDims(0);
			dims_data_dest[1] = carma_obj_handler_src->getDims(2);
			dims_data_dest[2] = carma_obj_handler_src->getDims(1);
			handle_dest->carma_object = new caObjD(context_handle,
					dims_data_dest);
			caObjD *carma_obj_handler_dest = (caObjD *) (handle_dest->carma_object);
			carma_obj_handler_dest->transpose(carma_obj_handler_src);
		}
	}
}

/*
 *                      _
 *  _ __ __ _ _ __   __| | ___  _ __ ___
 * | '__/ _` | '_ \ / _` |/ _ \| '_ ` _ \
 * | | | (_| | | | | (_| | (_) | | | | | |
 * |_|  \__,_|_| |_|\__,_|\___/|_| |_| |_|
 */

void Y_yoga_random(int argc) {
	/** @brief wrapper routine for yoga random method (uniform distribution)
	 *  @param[in] argc : command line arguments
	 *  can work as a (1) subroutine (return discarded) or (2) as a function
	 *    - first  : (1) the yoga_obj to be filled with random numbers
	 *                   / (2) type of object to be created
	 *    - second : (1) optional seed / (2) dimensions of object (Yorick convention)
	 *    - third  : (2) optional device id
	 *  in case (2) the destination is pushed on the stack as a yoga_obj
	 */
	int seed = 1234;
	if (yarg_subroutine()) {
		yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
		if (argc > 1)
			seed = ygets_i(argc-2);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle->device);
		if (handle->type == Y_FLOAT) {
			caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
			if (carma_obj_handler->is_rng_init() == false) {
				//carma_obj_handler->init_prng(handle->device);
				carma_obj_handler->init_prng_host(seed);
			}
			//carma_obj_handler->prng('U');
			carma_obj_handler->prng_host('U');
		} else if (handle->type == Y_DOUBLE) {
			caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
			if (carma_obj_handler->is_rng_init() == false) {
				//carma_obj_handler->init_prng(handle->device);
				carma_obj_handler->init_prng_host(seed);
			}
			//carma_obj_handler->prng('U');
			carma_obj_handler->prng_host('U');
		} else if (handle->type == Y_SCOMPLEX) {
			caObjC *carma_obj_handler = (caObjC *) (handle->carma_object);
			if (carma_obj_handler->is_rng_init() == false) {
				carma_obj_handler->init_prng(handle->device);
			}
			carma_obj_handler->prng('U');
		} else if (handle->type == Y_COMPLEX) {
			caObjZ *carma_obj_handler = (caObjZ *) (handle->carma_object);
			if (carma_obj_handler->is_rng_init() == false) {
				carma_obj_handler->init_prng(handle->device);
			}
			carma_obj_handler->prng('U');
		}
	} else {
		// called as a function : we need to create an object
		long ntot;
		long dims;

		int yType = yarg_typeid(argc - 1);

		if (yType == Y_STRING) {
			char *type_data = ygets_q(argc - 1);
			if (strcmp(type_data, "float") == 0) {
				yType = Y_FLOAT;
			} else if (strcmp(type_data, "double") == 0) {
				yType = Y_DOUBLE;
			} else if (strcmp(type_data, "scomplex") == 0) {
				yType = Y_SCOMPLEX;
			} else if (strcmp(type_data, "complex") == 0) {
				yType = Y_COMPLEX;
			}
		} else
			y_error("expecting a string for the type and a list of dimensions");
		long *dims_data = ygeta_l(argc - 2, &ntot, &dims);
		carma_context *context_handle = _getCurrentContext();
		int activeDevice = context_handle->get_activeDevice();
		int mdevice = activeDevice;
		if (argc > 2) {
			mdevice = ygets_i(argc-3);
		}

		activeDevice = context_handle->set_activeDevice(mdevice);

		yObj_struct *handle = (yObj_struct *) ypush_obj(&yObj,
				sizeof(yObj_struct));

		handle->device = mdevice;
		handle->type = yType;

		if (yType == Y_FLOAT) {
			handle->carma_object = new caObjS(context_handle, dims_data);
		} else if (yType == Y_DOUBLE) {
			handle->carma_object = new caObjD(context_handle, dims_data);
		} else if (yType == Y_SCOMPLEX) {
			handle->carma_object = new caObjC(context_handle, dims_data);
		} else if (yType == Y_COMPLEX) {
			handle->carma_object = new caObjZ(context_handle, dims_data);
		}

		if (handle->type == Y_FLOAT) {
			caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
			//carma_obj_handler->init_prng(handle->device);
			carma_obj_handler->init_prng_host(seed);
			//carma_obj_handler->prng('U');
			carma_obj_handler->prng_host('U');
		} else if (handle->type == Y_DOUBLE) {
			caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
			//carma_obj_handler->init_prng(handle->device);
			carma_obj_handler->init_prng_host(seed);
			//carma_obj_handler->prng('U');
			carma_obj_handler->prng_host('U');
		} else if (handle->type == Y_SCOMPLEX) {
			caObjC *carma_obj_handler = (caObjC *) (handle->carma_object);
			carma_obj_handler->init_prng(handle->device);
			carma_obj_handler->prng('U');
		} else if (handle->type == Y_COMPLEX) {
			caObjZ *carma_obj_handler = (caObjZ *) (handle->carma_object);
			carma_obj_handler->init_prng(handle->device);
			carma_obj_handler->prng('U');
		}
	}
	cutilSafeThreadSync();
}

void Y_yoga_random_n(int argc) {
	/** @brief wrapper routine for yoga random method (normal distribution)
	 *  @param[in] argc : command line arguments
	 *  can work as a (1) subroutine (return discarded) or (2) as a function
	 *    - first  : (1) the yoga_obj to be filled with random numbers
	 *                   / (2) type of object to be created
	 *    - second : (1) optional seed / (2) dimensions of object (Yorick convention)
	 *    - third  : (2) optional device id
	 *  in case (2) the destination is pushed on the stack as a yoga_obj
	 */
	int seed = 1234;
	if (yarg_subroutine()) {
		yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
		if (argc > 1)
			seed = ygets_i(argc-2);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle->device);
		if (handle->type == Y_FLOAT) {
			caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
			if (carma_obj_handler->is_rng_init() == false) {
				//carma_obj_handler->init_prng(handle->device);
				carma_obj_handler->init_prng_host(seed);
			}
			//carma_obj_handler->prng('N');
			carma_obj_handler->prng_host('N');
		} else if (handle->type == Y_DOUBLE) {
			caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
			if (carma_obj_handler->is_rng_init() == false) {
				//carma_obj_handler->init_prng(handle->device);
				carma_obj_handler->init_prng_host(seed);
			}
			//carma_obj_handler->prng('N');
			carma_obj_handler->prng_host('N');
		} else if (handle->type == Y_SCOMPLEX) {
			caObjC *carma_obj_handler = (caObjC *) (handle->carma_object);
			if (carma_obj_handler->is_rng_init() == false) {
				carma_obj_handler->init_prng(handle->device);
			}
			carma_obj_handler->prng('N');
		} else if (handle->type == Y_COMPLEX) {
			caObjZ *carma_obj_handler = (caObjZ *) (handle->carma_object);
			if (carma_obj_handler->is_rng_init() == false) {
				carma_obj_handler->init_prng(handle->device);
			}
			carma_obj_handler->prng('N');
		}
	} else {
		// called as a function : we need to create an object
		long ntot;
		long dims;

		int yType = yarg_typeid(argc - 1);
		if (yType == Y_STRING) {
			char *type_data = ygets_q(argc - 1);
			if (strcmp(type_data, "float") == 0) {
				yType = Y_FLOAT;
			} else if (strcmp(type_data, "double") == 0) {
				yType = Y_DOUBLE;
			} else if (strcmp(type_data, "scomplex") == 0) {
				yType = Y_SCOMPLEX;
			} else if (strcmp(type_data, "complex") == 0) {
				yType = Y_COMPLEX;
			}
		}

		long *dims_data = ygeta_l(argc - 2, &ntot, &dims);

		carma_context *context_handle = _getCurrentContext();
		int activeDevice = context_handle->get_activeDevice();

		int mdevice = activeDevice;
		if (argc > 2) {
			mdevice = ygets_i(argc-3);
		}
		activeDevice = context_handle->set_activeDevice(mdevice);

		yObj_struct *handle = (yObj_struct *) ypush_obj(&yObj,
				sizeof(yObj_struct));
		handle->device = mdevice;
		handle->type = yType;
		if (yType == Y_FLOAT) {
			handle->carma_object = new caObjS(context_handle, dims_data);
		} else if (yType == Y_DOUBLE) {
			handle->carma_object = new caObjD(context_handle, dims_data);
		} else if (yType == Y_SCOMPLEX) {
			handle->carma_object = new caObjC(context_handle, dims_data);
		} else if (yType == Y_COMPLEX) {
			handle->carma_object = new caObjZ(context_handle, dims_data);
		}

		if (handle->type == Y_FLOAT) {
			caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
			//carma_obj_handler->init_prng(handle->device);
			carma_obj_handler->init_prng_host(seed);
			//carma_obj_handler->prng('N');
			carma_obj_handler->prng_host('N');
		} else if (handle->type == Y_DOUBLE) {
			caObjD *carma_obj_handler = (caObjD *) (handle->carma_object);
			//carma_obj_handler->init_prng(handle->device);
			carma_obj_handler->init_prng_host(seed);
			//carma_obj_handler->prng('N');
			carma_obj_handler->prng_host('N');
		} else if (handle->type == Y_SCOMPLEX) {
			caObjC *carma_obj_handler = (caObjC *) (handle->carma_object);
			carma_obj_handler->init_prng(handle->device);
			carma_obj_handler->prng('N');
		} else if (handle->type == Y_COMPLEX) {
			caObjZ *carma_obj_handler = (caObjZ *) (handle->carma_object);
			carma_obj_handler->init_prng(handle->device);
			carma_obj_handler->prng('N');
		}
	}
	cutilSafeThreadSync();
}

/*
 *    __  __ _
 *   / _|/ _| |_
 *  | |_| |_| __|
 *  |  _|  _| |_
 *  |_| |_|  \__|
 */

void Y_yoga_fft(int argc) {
	/** @brief wrapper routine for yoga fft method
	 *  @param[in] argc : command line arguments
	 *  can work as a (1) subroutine (return discarded) or (2) as a function
	 *    - first  : (1) / (2) the yoga_obj to tranform
	 *    - second : (1) optional yoga_obj to contain transform (out of place case)
	 *                   / (2) optional direction
	 *    - third  : (1) optional direction
	 *  in case (2) the transform is out of place and result is pushed on the stack as a yoga_obj
	 */
	if (yarg_subroutine()) {
		yObj_struct *handle_src = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		if (argc > 1) {
			// out of place
			yObj_struct *handle_dest = (yObj_struct *) yget_obj(argc - 2,
					&yObj);
			if (handle_dest->device != handle_src->device)
				y_error("transpose only on the same device");
			carma_context *context_handle = _getCurrentContext();
			context_handle->set_activeDevice(handle_dest->device);
			int dir;
			if (argc > 3)
				dir = ygets_i(argc-3);
			else
				dir = 1;

			if (handle_src->type == Y_FLOAT) {
				if (handle_dest->type != Y_SCOMPLEX)
					throw "destination has wrong type";
				else {
					caObjS *carma_obj_handler_src =
							(caObjS *) (handle_src->carma_object);
					caObjC *carma_obj_handler_dest =
							(caObjC *) (handle_dest->carma_object);

					if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<float, cuFloatComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<float, cuFloatComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					  /*
						carma_initfft<float, cuFloatComplex>(
								carma_obj_handler_src->getDims(),
								*carma_obj_handler_src->getPlan(),
								carma_select_plan<float, cuFloatComplex>());
					  */
					}
					carma_fft(carma_obj_handler_src->getData(),
							carma_obj_handler_dest->getData(), dir,
							*carma_obj_handler_src->getPlan());
				}
			}
			if (handle_src->type == Y_SCOMPLEX) {
				if (handle_dest->type == Y_SCOMPLEX) {
					caObjC *carma_obj_handler_src =
							(caObjC *) (handle_src->carma_object);
					caObjC *carma_obj_handler_dest =
							(caObjC *) (handle_dest->carma_object);
					if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<cuFloatComplex, cuFloatComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<cuFloatComplex, cuFloatComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
						carma_initfft<cuFloatComplex, cuFloatComplex>(
								carma_obj_handler_src->getDims(),
								*carma_obj_handler_src->getPlan(),
								carma_select_plan<cuFloatComplex, cuFloatComplex>());
					    */
					}
					carma_fft(carma_obj_handler_src->getData(),
							carma_obj_handler_dest->getData(), dir,
							*carma_obj_handler_src->getPlan());
				} else if (handle_dest->type == Y_FLOAT) {
					caObjC *carma_obj_handler_src =
							(caObjC *) (handle_src->carma_object);
					caObjS *carma_obj_handler_dest =
							(caObjS *) (handle_dest->carma_object);
					if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<cuFloatComplex, float>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<cuFloatComplex, float>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
						carma_initfft<cuFloatComplex, float>(
								carma_obj_handler_src->getDims(),
								*carma_obj_handler_src->getPlan(),
								carma_select_plan<cuFloatComplex, float>());
					    */
					}
					carma_fft(carma_obj_handler_src->getData(),
							carma_obj_handler_dest->getData(), dir,
							*carma_obj_handler_src->getPlan());
				} else
					throw "destination has wrong type";
			} else if (handle_src->type == Y_DOUBLE) {
				if (handle_dest->type == Y_COMPLEX) {
					caObjD *carma_obj_handler_src =
							(caObjD *) (handle_src->carma_object);
					caObjZ *carma_obj_handler_dest =
							(caObjZ *) (handle_dest->carma_object);
					if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<double, cuDoubleComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<double, cuDoubleComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
						carma_initfft<double, cuDoubleComplex>(
								carma_obj_handler_src->getDims(),
								*carma_obj_handler_src->getPlan(),
								carma_select_plan<double, cuDoubleComplex>());
					    */
					}
					carma_fft(carma_obj_handler_src->getData(),
							carma_obj_handler_dest->getData(), dir,
							*carma_obj_handler_src->getPlan());
				} else
					throw "destination has wrong type";
			} else if (handle_src->type == Y_COMPLEX) {
				if (handle_dest->type == Y_COMPLEX) {
					caObjZ *carma_obj_handler_src =
							(caObjZ *) (handle_src->carma_object);
					caObjZ *carma_obj_handler_dest =
							(caObjZ *) (handle_dest->carma_object);
					if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<cuDoubleComplex, cuDoubleComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<cuDoubleComplex, cuDoubleComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
						carma_initfft<cuDoubleComplex, cuDoubleComplex>(
								carma_obj_handler_src->getDims(),
								*carma_obj_handler_src->getPlan(),
								carma_select_plan<cuDoubleComplex,
										cuDoubleComplex>());
					    */
					}
					carma_fft(carma_obj_handler_src->getData(),
							carma_obj_handler_dest->getData(), dir,
							*carma_obj_handler_src->getPlan());
				} else if (handle_dest->type == Y_DOUBLE) {
					caObjZ *carma_obj_handler_src =
							(caObjZ *) (handle_src->carma_object);
					caObjD *carma_obj_handler_dest =
							(caObjD *) (handle_dest->carma_object);
					if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<cuDoubleComplex, double>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<cuDoubleComplex, double>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
						carma_initfft<cuDoubleComplex, double>(
								carma_obj_handler_src->getDims(),
								*carma_obj_handler_src->getPlan(),
								carma_select_plan<cuDoubleComplex, double>());
					    */
					}
					carma_fft(carma_obj_handler_src->getData(),
							carma_obj_handler_dest->getData(), dir,
							*carma_obj_handler_src->getPlan());
				} else
					throw "destination has wrong type";
			} else
				throw "wrong type for fft";
		} else {
			int dir;
			carma_context *context_handle = _getCurrentContext();
			context_handle->set_activeDevice(handle_src->device);
			if (argc > 1)
				dir = ygets_i(argc-2);
			else
				dir = 1;

			if (handle_src->type == Y_SCOMPLEX) {
				caObjC *carma_obj_handler_src =
						(caObjC *) (handle_src->carma_object);
				if (*carma_obj_handler_src->getPlan() == 0L) {				  
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<cuFloatComplex, cuFloatComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<cuFloatComplex, cuFloatComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
					carma_initfft<cuFloatComplex, cuFloatComplex>(
							carma_obj_handler_src->getDims(),
							*plan,
							carma_select_plan<cuFloatComplex, cuFloatComplex>());
					    */
				}
				carma_fft(carma_obj_handler_src->getData(),
						carma_obj_handler_src->getData(), dir,
						*carma_obj_handler_src->getPlan());
			} else if (handle_src->type == Y_COMPLEX) {
				caObjD *carma_obj_handler_src =
						(caObjD *) (handle_src->carma_object);
				if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<cuDoubleComplex, cuDoubleComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<cuDoubleComplex, cuDoubleComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
					carma_initfft<cuDoubleComplex, cuDoubleComplex>(
							carma_obj_handler_src->getDims(),
							*carma_obj_handler_src->getPlan(),
							carma_select_plan<cuDoubleComplex, cuDoubleComplex>());
					    */
				}
				carma_fft(carma_obj_handler_src->getData(),
						carma_obj_handler_src->getData(), dir,
						*carma_obj_handler_src->getPlan());
			} else
				throw "wrong type for fft";
		}
	} else {
		yObj_struct *handle_src = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_src->device);
		// called as a function : we need to create an object to push on the stack
		// check if src plan exists
		// create if not and execute
		yObj_struct *handle = (yObj_struct *) ypush_obj(&yObj,
				sizeof(yObj_struct));
		handle->device = handle_src->device;
		int dir;
		if (argc > 2)
			dir = ygets_i(argc-2);
		else
			dir = 1;

		if ((handle_src->type == Y_FLOAT) || (handle_src->type == Y_SCOMPLEX)) {
			if (handle_src->type == Y_FLOAT) {
				caObjS *carma_obj_handler_src =
						(caObjS *) (handle_src->carma_object);
				if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<float, cuFloatComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<float, cuFloatComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
					carma_initfft<float, cuFloatComplex>(
							carma_obj_handler_src->getDims(),
							*carma_obj_handler_src->getPlan(),
							carma_select_plan<float, cuFloatComplex>());
					    */
				}
				handle->carma_object = new caObjC(context_handle,
						carma_obj_handler_src->getDims());
				handle->type = Y_SCOMPLEX;
				caObjC *carma_obj_handler = (caObjC *) (handle->carma_object);
				carma_fft(carma_obj_handler_src->getData(),
						carma_obj_handler->getData(), dir,
						*carma_obj_handler_src->getPlan());
			} else {
				caObjC *carma_obj_handler_src =
						(caObjC *) (handle_src->carma_object);
				if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<cuFloatComplex, cuFloatComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<cuFloatComplex, cuFloatComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
					carma_initfft<cuFloatComplex, cuFloatComplex>(
							carma_obj_handler_src->getDims(),
							*carma_obj_handler_src->getPlan(),
							carma_select_plan<cuFloatComplex, cuFloatComplex>());
					    */
				}
				handle->carma_object = new caObjC(context_handle,
						carma_obj_handler_src->getDims());
				handle->type = Y_SCOMPLEX;
				caObjC *carma_obj_handler = (caObjC *) (handle->carma_object);
				carma_fft(carma_obj_handler_src->getData(),
						carma_obj_handler->getData(), dir,
						*carma_obj_handler_src->getPlan());
			}
		} else if ((handle_src->type == Y_DOUBLE)
				|| (handle_src->type == Y_COMPLEX)) {
			if (handle_src->type == Y_DOUBLE) {
				caObjD *carma_obj_handler_src =
						(caObjD *) (handle_src->carma_object);
				if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<double, cuDoubleComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<double, cuDoubleComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
					carma_initfft<double, cuDoubleComplex>(
							carma_obj_handler_src->getDims(),
							*carma_obj_handler_src->getPlan(),
							carma_select_plan<double, cuDoubleComplex>());
					    */
				}
				handle->carma_object = new caObjZ(context_handle,
						carma_obj_handler_src->getDims());
				handle->type = Y_COMPLEX;
				caObjZ *carma_obj_handler = (caObjZ *) (handle->carma_object);
				carma_fft(carma_obj_handler_src->getData(),
						carma_obj_handler->getData(), dir,
						*carma_obj_handler_src->getPlan());
			} else {
				caObjZ *carma_obj_handler_src =
						(caObjZ *) (handle_src->carma_object);
				if (*carma_obj_handler_src->getPlan() == 0L) {
					  cufftHandle *plan=carma_obj_handler_src->getPlan();///< FFT plan
					  if(carma_obj_handler_src->getDims()[0]==2)
					    /* Create a 2D FFT plan. */ 
					    cufftSafeCall( cufftPlan2d(plan, carma_obj_handler_src->getDims()[1], 
								       carma_obj_handler_src->getDims()[2], 
								       carma_select_plan<cuDoubleComplex, cuDoubleComplex>()));
					  else
					    /* Create a 3D FFT plan. */ {
					    int mdims[2];
					    mdims[0] = (int)(carma_obj_handler_src->getDims()[1]);
					    mdims[1] = (int)(carma_obj_handler_src->getDims()[2]);
					    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, 
									 carma_select_plan<cuDoubleComplex, cuDoubleComplex>(),
									 (int)(carma_obj_handler_src->getDims()[3])));
					  }
					    /*
					carma_initfft<cuDoubleComplex, cuDoubleComplex>(
							carma_obj_handler_src->getDims(),
							*carma_obj_handler_src->getPlan(),
							carma_select_plan<cuDoubleComplex, cuDoubleComplex>());
					    */
				}
				handle->carma_object = new caObjZ(context_handle,
						carma_obj_handler_src->getDims());
				handle->type = Y_COMPLEX;
				caObjZ *carma_obj_handler = (caObjZ *) (handle->carma_object);
				carma_fft(carma_obj_handler_src->getData(),
						carma_obj_handler->getData(), dir,
						*carma_obj_handler_src->getPlan());
			}
		}
	}
	cutilSafeThreadSync();
}

/*
 *  ___  ___ __ _ _ __
 * / __|/ __/ _` | '_ \
 * \__ \ (_| (_| | | | |
 * |___/\___\__,_|_| |_|

 */

void Y_yoga_max(int argc)
/** @brief wrapper routine for yoga max method (see cudpp library)
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  result is pushed on stack
 *  only single precision int and float are supported
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(handle->device);
	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		carma_obj_handler->config_scan((char *) ("max"), 1, 1);
		ypush_double((double) carma_obj_handler->scan());
	}
	if (handle->type == Y_INT) {
		caObjI *carma_obj_handler = (caObjI *) (handle->carma_object);
		carma_obj_handler->config_scan((char *) ("max"), 1, 1);
		ypush_int(carma_obj_handler->scan());
	}
	if (handle->type == Y_DOUBLE) {
		throw "Double not supported";
	}
}

void Y_yoga_min(int argc)
/** @brief wrapper routine for yoga min method (see cudpp library)
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  result is pushed on stack
 *  only single precision int and float are supported
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(handle->device);
	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		carma_obj_handler->config_scan((char *) ("min"), 1, 1);
		ypush_double((double) carma_obj_handler->scan());
	}
	if (handle->type == Y_INT) {
		caObjI *carma_obj_handler = (caObjI *) (handle->carma_object);
		carma_obj_handler->config_scan((char *) ("min"), 1, 1);
		ypush_int(carma_obj_handler->scan());
	}
	if (handle->type == Y_DOUBLE) {
		throw "Double not supported";
	}
}

void Y_yoga_mult(int argc)
/** @brief wrapper routine for yoga mult method (see cudpp library)
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  result is pushed on stack
 *  only single precision int and float are supported
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(handle->device);
	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		carma_obj_handler->config_scan((char *) ("mult"), 1, 1);
		ypush_double((double) carma_obj_handler->scan());
	}
	if (handle->type == Y_INT) {
		caObjI *carma_obj_handler = (caObjI *) (handle->carma_object);
		carma_obj_handler->config_scan((char *) ("mult"), 1, 1);
		ypush_int(carma_obj_handler->scan());
	}
	if (handle->type == Y_DOUBLE) {
		throw "Double not supported";
	}
}

void Y_yoga_add(int argc)
/** @brief wrapper routine for yoga add method (see cudpp library)
 *  @param[in] argc : command line argument (yoga_obj expected)
 *  result is pushed on stack
 *  only single precision int and float are supported
 */
{
	yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDevice(handle->device);
	if (handle->type == Y_FLOAT) {
		caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
		carma_obj_handler->config_scan((char *) ("add"), 1, 1);
		ypush_double((double) carma_obj_handler->scan());
	}
	if (handle->type == Y_INT) {
		caObjI *carma_obj_handler = (caObjI *) (handle->carma_object);
		carma_obj_handler->config_scan((char *) ("add"), 1, 1);
		ypush_int(carma_obj_handler->scan());
	}
	if (handle->type == Y_DOUBLE) {
		throw "Double not supported";
	}
}

/*
 *                 _
 *  ___  ___  _ __| |_
 * / __|/ _ \| '__| __|
 * \__ \ (_) | |  | |_
 * |___/\___/|_|   \__|
 */

void Y_yoga_sort(int argc)
/** @brief wrapper routine for yoga sort method
 *  @param[in] argc : command line arguments
 *  can work as a (1) subroutine (return discarded) or (2) as a function
 *    - first  : (1) / (2) the yoga_obj of keys to be sorted
 *    - second : (1) optional object to handle values
 *  in case (2) the sorted array is pushed on the stack as a yoga_obj
 */
{
	if (yarg_subroutine()) {
		bool keysOnly;
		void *val_handle = NULL;
		void *carma_val_handler = NULL;
		yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDeviceForCpy(handle->device);
		if (handle->type == Y_FLOAT) {
			if (argc > 1) {
				val_handle = (yObj_struct *) yget_obj(argc - 2, &yObj);
				if (handle->device != ((yObj_struct *) val_handle)->device)
					y_error("keys & val for sort should be on same device");
				keysOnly = false;
			} else
				keysOnly = true;
			caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
			if (!keysOnly) {
				carma_val_handler =
						(caObjUI *) (((yObj_struct *) val_handle)->carma_object);
				carma_obj_handler->device2deviceInd(
						((caObjUI *) carma_val_handler)->getData());
			}
			carma_obj_handler->sort_init(keysOnly);
			carma_obj_handler->sort();
			if (!keysOnly) {
				cudaMemcpy(((caObjUI *) carma_val_handler)->getData(),
						((caObjS *) carma_obj_handler)->getValues(),
						((caObjS *) carma_obj_handler)->getNbElem()
								* sizeof(unsigned int),
						cudaMemcpyDeviceToDevice);
			}
		} else if (handle->type == Y_DOUBLE) {
			throw "Double not supported";
		}
	} else {
		// called as a function : we need to create an object
		yObj_struct *handle = (yObj_struct *) yget_obj(argc - 1, &yObj);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle->device);
		yObj_struct *val_handle = (yObj_struct *) ypush_obj(&yObj,
				sizeof(yObj_struct));
		((yObj_struct *) val_handle)->device = handle->device;
		void *carma_val_handler;
		val_handle->type = Y_SHORT;
		if (handle->type == Y_FLOAT) {
			caObjS *carma_obj_handler = (caObjS *) (handle->carma_object);
			val_handle->carma_object = new caObjUI(context_handle,
					carma_obj_handler->getDims());
			carma_val_handler =
					(caObjUI *) (((yObj_struct *) val_handle)->carma_object);
			carma_obj_handler->sort_init(false);
			carma_obj_handler->sort();
			cudaMemcpy(((caObjUI *) carma_val_handler)->getData(),
					((caObjS *) carma_obj_handler)->getValues(),
					((caObjS *) carma_obj_handler)->getNbElem()
							* sizeof(unsigned int), cudaMemcpyDeviceToDevice);
		} else if (handle->type == Y_DOUBLE) {
			throw "Double not supported";
		}
	}
}

/*
 *                                       _
 *   ___ ___  _ __ ___  _ __   __ _  ___| |_
 *  / __/ _ \| '_ ` _ \| '_ \ / _` |/ __| __|
 * | (_| (_) | | | | | | |_) | (_| | (__| |_
 *  \___\___/|_| |_| |_| .__/ \__,_|\___|\__|
 *                     |_|
 */

void Y_yoga_compact(int argc)
/** @brief wrapper routine for yoga compact method (see cudpp library)
 *  @param[in] argc : command line arguments
 *  can work as a (1) subroutine (return discarded) or (2) as a function
 *    - first  : (1) / (2) the yoga_obj to be compacted
 *    - second : (1) the compacted yoga_obj
 *                   / (2) valid indexes can be either a standard array or a yoga_obj
 *    - third  : (1) optional valid indexes can be either a standard array or a yoga_obj
 *  in case (2) the compacted array is pushed on the stack as a yoga_obj
 *  only single precision floating point is supported
 *  valid indexes can have the dimension of dest or src
 *  if it has the dim of dest, it is assumed as a list of valid indexes
 *  if it has the dim of src it is assumed as an array of 0 & 1 with 1 for valids
 *  in case (2), valid indexes must have the dim of dest
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle_src = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_dest = (yObj_struct *) yget_obj(argc - 2,
				&yObj);

		if (handle_src->device != handle_dest->device)
			y_error("src & dest for compact should be on same device");

		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDeviceForCpy(handle_src->device);

		int yType;
		unsigned int *isValid;
		long int nb_valids;

		if (argc > 2) {
			yType = yarg_typeid(argc - 3);
			if (handle_dest->type == Y_FLOAT) {
				caObjS *carma_src_handler = (caObjS *) (handle_src->carma_object);
				caObjS *carma_dest_handler = (caObjS *) (handle_dest->carma_object);

				// check if sort_on on src. if not initialize values (needed)
				unsigned int *values = carma_src_handler->getValues();
				// check if sort_on on src. if not initialize values (needed)
				if (values == NULL) {
					cutilSafeCall(
							cudaMalloc((void **)&values,carma_src_handler->getNbElem() * sizeof(unsigned int)));
				}
				cutilSafeCall(
						cudaMemset(values, 0, sizeof(unsigned int)*carma_src_handler->getNbElem()));

				if (yType == Y_OPAQUE) {
					// the user gave a carma_obj as an array of valids
					yObj_struct *handle_valids = (yObj_struct *) yget_obj(
							argc - 3, &yObj);
					if (handle_valids->device != handle_dest->device)
						y_error(
								"valids & dest for compact should be on same device");
					caObjUI *carma_valids_handler =
							(caObjUI *) (handle_valids->carma_object);
					nb_valids = (long) carma_valids_handler->getNbElem();
					if (nb_valids != carma_src_handler->getNbElem()) {
						fillvalues(carma_src_handler->getValues(),
								carma_valids_handler->getData(),
								carma_valids_handler->getNbElem());
					} else {
						cutilSafeCall(
								cudaMemcpy(carma_src_handler->getValues(), carma_valids_handler->getData(), carma_src_handler->getNbElem() * sizeof(unsigned int), cudaMemcpyDeviceToDevice));
					}
				} else {
					// the user gave a yorick array as an array of valids
					long dims[Y_DIMSIZE];
					isValid = (unsigned int *) ygeta_i(argc - 3, &nb_valids,
							dims);
					if (nb_valids != carma_src_handler->getNbElem()) {
						unsigned int *d_valids;
						cutilSafeCall(
								cudaMalloc((void **)&d_valids,nb_valids * sizeof(unsigned int)));
						cudaMemcpy(d_valids, isValid,
								nb_valids * sizeof(unsigned int),
								cudaMemcpyHostToDevice);
						fillvalues(carma_src_handler->getValues(), d_valids,
								nb_valids);
						cudaFree(d_valids);
					} else {
						cudaMemcpy(carma_src_handler->getValues(), isValid,
								carma_src_handler->getNbElem()
										* sizeof(unsigned int),
								cudaMemcpyDeviceToDevice);
					}
				}
				carma_src_handler->compact_init();
				carma_src_handler->compact(carma_dest_handler);
			} else if (handle_src->type == Y_DOUBLE) {
				throw "Double not supported";
			}
		} else {
			//we assume everything is already set we just compact
			if (handle_dest->type == Y_FLOAT) {
				caObjS *carma_src_handler = (caObjS *) (handle_src->carma_object);
				caObjS *carma_dest_handler = (caObjS *) (handle_dest->carma_object);
				carma_src_handler->compact(carma_dest_handler);
			} else if (handle_src->type == Y_DOUBLE) {
				throw "Double not supported";
			}
		}
	} else {
		// called as a function : we need to create an object
		yObj_struct *handle_src = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_src->device);
		int yType = yarg_typeid(argc - 2);
		unsigned int *isValid;
		long int nb_valids;
		if (handle_src->type == Y_FLOAT) {
			caObjS *carma_src_handler = (caObjS *) (handle_src->carma_object);

			unsigned int *values = carma_src_handler->getValues();
			// check if sort_on on src. if not initialize values (needed)
			if (values == NULL) {
				cutilSafeCall(
						cudaMalloc((void **)&values,carma_src_handler->getNbElem() * sizeof(unsigned int)));
			}
			cutilSafeCall(
					cudaMemset(values, 0, sizeof(unsigned int)*carma_src_handler->getNbElem()));

			if (yType == Y_OPAQUE) {
				// the user gave a carma_obj as an array of valids
				yObj_struct *handle_valids = (yObj_struct *) yget_obj(argc - 2,
						&yObj);
				caObjUI *carma_valids_handler =
						(caObjUI *) (handle_valids->carma_object);
				if (handle_valids->device != handle_src->device)
					y_error(
							"valids & src for compact should be on same device");
				nb_valids = carma_valids_handler->getNbElem();
				if (nb_valids != carma_src_handler->getNbElem()) {
					fillvalues(carma_src_handler->getValues(),
							carma_valids_handler->getData(),
							carma_valids_handler->getNbElem());
				} else {
					throw "operation not supported";
				}
			} else {
				// the user gave a yorick array as an array of valids
				long dims[Y_DIMSIZE];
				isValid = (unsigned int *) ygeta_i(argc - 2, &nb_valids, dims);
				if (nb_valids != carma_src_handler->getNbElem()) {
					unsigned int *d_valids;
					cutilSafeCall(
							cudaMalloc((void **)&d_valids,nb_valids * sizeof(unsigned int)));
					cudaMemcpy(d_valids, isValid,
							nb_valids * sizeof(unsigned int),
							cudaMemcpyHostToDevice);
					fillvalues(carma_src_handler->getValues(), d_valids,
							nb_valids);
					cudaFree(d_valids);
				} else {
					throw "operation not supported";
				}
			}
			long dims_data[2];
			dims_data[0] = 1;
			dims_data[1] = (long) nb_valids;
			yObj_struct *handle_dest = (yObj_struct *) ypush_obj(&yObj,
					sizeof(yObj_struct));
			handle_dest->device = handle_src->device;
			handle_dest->type = Y_FLOAT;
			handle_dest->carma_object = new caObjS(context_handle, dims_data);
			caObjS * carma_dest_handler = (caObjS *) (handle_dest->carma_object);
			carma_src_handler->compact_init();
			carma_src_handler->compact(carma_dest_handler);
		} else if (handle_src->type == Y_DOUBLE) {
			throw "Double not supported";
		}
	}
}

/*
 *   __  __ _
 *  / _|/ _| |_     ___ ___  _ ____   __
 * | |_| |_| __|   / __/ _ \| '_ \ \ / /
 * |  _|  _| |_   | (_| (_) | | | \ V /
 * |_| |_|  \__|___\___\___/|_| |_|\_/
 *            |_____|
 */

void Y_yoga_fftconv(int argc)
/** @brief wrapper routine for fft convolution
 *  @param[in] argc : command line arguments
 *  can only be called as a subroutine (return discarded)
 *    - first  : the yoga_obj containing the result
 *    - second : the yoga_obj containing the image
 *    - third  : the yoga_obj containing the kernel
 *    - fourth : the yoga_obj containing the padded_data
 *    - fifth  : the yoga_obj containing the padded_kernel
 *  5 yoga_obj must be initialized : image, kernel, padded_data, padded_kernel and result
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle_res = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_im = (yObj_struct *) yget_obj(argc - 2, &yObj);
		yObj_struct *handle_ker = (yObj_struct *) yget_obj(argc - 3,
				&yObj);
		yObj_struct *handle_pdata = (yObj_struct *) yget_obj(argc - 4,
				&yObj);
		yObj_struct *handle_pspec = (yObj_struct *) yget_obj(argc - 5,
				&yObj);
		if ((handle_res->device != handle_im->device)
				|| (handle_res->device != handle_ker->device)
				|| (handle_res->device != handle_pdata->device)
				|| (handle_res->device != handle_pspec->device))
			y_error("fftconv only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_res->device);
		int kernelY = ygets_i(argc-6);
		int kernelX = ygets_i(argc-7);

		if ((handle_res->type == Y_FLOAT) && (handle_im->type == Y_FLOAT)
				&& (handle_ker->type == Y_FLOAT)
				&& (handle_pdata->type == Y_FLOAT)
				&& (handle_pspec->type == Y_SCOMPLEX)) {
			caObjS *carma_res_handler = (caObjS *) (handle_res->carma_object);
			caObjS *carma_im_handler = (caObjS *) (handle_im->carma_object);
			caObjS *carma_ker_handler = (caObjS *) (handle_ker->carma_object);
			caObjS *carma_pdata_handler = (caObjS *) (handle_pdata->carma_object);
			caObjC *carma_pspec_handler = (caObjC *) (handle_pspec->carma_object);
			carma_initfftconv(carma_im_handler, carma_ker_handler,
					carma_pdata_handler, carma_pspec_handler, kernelY, kernelX);

			carma_fftconv(carma_res_handler, carma_pdata_handler,
					carma_pspec_handler, kernelY, kernelX);
		} else
			throw "args must be float float float float scomplex";
	}
}

void Y_yoga_fftconv_init(int argc)
/** @brief wrapper routine for fft convolution initialization
 *  @param[in] argc : command line arguments
 *  can only be called as a function
 *    - first  : the yoga_obj containing the image
 *    - second : the yoga_obj containing the result
 *    - third  : type of workspace to initialize (real or complex, i.e. padded or not)
 *  push the corresponding result object on the stack
 */
{
	if (yarg_subroutine())
		throw "can only be called as a function";
	else {
		yObj_struct *handle_im = (yObj_struct *) yget_obj(argc - 1, &yObj);
		yObj_struct *handle_ker = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		if (handle_im->device != handle_ker->device)
			y_error("fftconv only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_im->device);
		char *type_data = ygets_q(argc - 3);
		caObjS *carma_im_handler = (caObjS *) (handle_im->carma_object);
		caObjS *carma_ker_handler = (caObjS *) (handle_ker->carma_object);
		yObj_struct *handle_dest = (yObj_struct *) ypush_obj(&yObj,
				sizeof(yObj_struct));
		handle_dest->device = handle_im->device;

		const int dataH = carma_im_handler->getDims(1);
		const int dataW = carma_im_handler->getDims(2);
		const int kernelH = carma_ker_handler->getDims(1);
		const int kernelW = carma_ker_handler->getDims(2);

		const int fftH = snapTransformSize(dataH + kernelH - 1);
		const int fftW = snapTransformSize(dataW + kernelW - 1);

		long dims_data[carma_im_handler->getDims(0) + 1];
		dims_data[0] = carma_im_handler->getDims(0);

		if (strcmp(type_data, "real") == 0) {
			dims_data[1] = fftH;
			dims_data[2] = fftW;
			if (carma_im_handler->getDims(0) > 2)
				dims_data[3] = carma_im_handler->getDims(3);
			handle_dest->type = Y_FLOAT;
			handle_dest->device = handle_im->device;
			handle_dest->carma_object = new caObjS(context_handle, dims_data);
		} else if (strcmp(type_data, "complex") == 0) {
			dims_data[1] = fftH;
			dims_data[2] = fftW / 2 + 1;
			if (carma_im_handler->getDims(0) > 2)
				dims_data[3] = carma_im_handler->getDims(3);
			handle_dest->type = Y_SCOMPLEX;
			handle_dest->device = handle_im->device;
			handle_dest->carma_object = new caObjC(context_handle, dims_data);
		} else
			throw "fftconv ws type not supported";
	}
}

/*
 *                              _               _
 *  _   _  ___   __ _  __ _    | |__   ___  ___| |_
 * | | | |/ _ \ / _` |/ _` |   | '_ \ / _ \/ __| __|
 * | |_| | (_) | (_| | (_| |   | | | | (_) \__ \ |_
 *  \__, |\___/ \__, |\__,_|___|_| |_|\___/|___/\__|
 *  |___/       |___/     |_____|
 */

void yHostObj_free(void *obj)
/** @brief yHostObj_struct destructor.
 *  @param obj : yHostObj_struct to free
 */
{
	yHostObj_struct *handler = (yHostObj_struct *) obj;
	try {
		if (handler->type == Y_FLOAT) {
			carma_host_obj<float> *carma_host_obj_handler =
					(carma_host_obj<float> *) (handler->carma_host_object);
			delete carma_host_obj_handler;
		} else if (handler->type == Y_DOUBLE) {
			carma_host_obj<double> *carma_host_obj_handler =
					(carma_host_obj<double> *) (handler->carma_host_object);
			delete carma_host_obj_handler;
		} else if (handler->type == Y_INT) {
			carma_host_obj<int> *carma_host_obj_handler =
					(carma_host_obj<int> *) (handler->carma_host_object);
			delete carma_host_obj_handler;
		} else if (handler->type == Y_SHORT) {
			carma_host_obj<unsigned int> *carma_host_obj_handler = (carma_host_obj<
					unsigned int> *) (handler->carma_host_object);
			delete carma_host_obj_handler;
		} else if (handler->type == Y_COMPLEX) {
			carma_host_obj<cuFloatComplex> *carma_host_obj_handler =
					(carma_host_obj<cuFloatComplex> *) (handler->carma_host_object);
			delete carma_host_obj_handler;
		} else if (handler->type == Y_SCOMPLEX) {
			carma_host_obj<cuDoubleComplex> *carma_host_obj_handler =
					(carma_host_obj<cuDoubleComplex> *) (handler->carma_host_object);
			delete carma_host_obj_handler;
		} else
			throw "Type unknown";
		handler->type = Y_VOID;
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	}
}

void yHostObj_print(void *obj)
/** @brief yHostObj_struct printer.
 *  @param[in] obj : yHostObj_struct to print
 */
{
	ostringstream mystr;
	yHostObj_struct *handler = (yHostObj_struct *) obj;
	try {
		if (handler->type == Y_FLOAT) {
			carma_host_obj<float> *carma_host_obj_handler =
					(carma_host_obj<float> *) (handler->carma_host_object);
			mystr << "Carma Host Object: " << endl;
			mystr << "  type: float";
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			mystr << "  Allocation: " << carma_host_obj_handler->getMetAlloc()
					<< endl;
			mystr << "  nb streams: " << carma_host_obj_handler->get_nbStreams();
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			long *dims = carma_host_obj_handler->getDims();
			mystr << "  Dims: " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_DOUBLE) {
			carma_host_obj<double> *carma_host_obj_handler =
					(carma_host_obj<double> *) (handler->carma_host_object);
			mystr << "Carma Host Object : " << endl;
			mystr << "Type : double";
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			mystr << "Allocation : " << carma_host_obj_handler->getMetAlloc()
					<< endl;
			mystr << "Nb streams : " << carma_host_obj_handler->get_nbStreams();
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			long *dims = carma_host_obj_handler->getDims();
			mystr << "Dims : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_INT) {
			carma_host_obj<int> *carma_host_obj_handler =
					(carma_host_obj<int> *) (handler->carma_host_object);
			mystr << "Carma Host Object : " << endl;
			mystr << "Type : int";
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			mystr << "Allocation : " << carma_host_obj_handler->getMetAlloc()
					<< endl;
			mystr << "Nb streams : " << carma_host_obj_handler->get_nbStreams();
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			long *dims = carma_host_obj_handler->getDims();
			mystr << "Dims : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_SHORT) {
			carma_host_obj<unsigned int> *carma_host_obj_handler = (carma_host_obj<
					unsigned int> *) (handler->carma_host_object);
			mystr << "Carma Host Object : " << endl;
			mystr << "Type : unsigned int";
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			mystr << "Allocation : " << carma_host_obj_handler->getMetAlloc()
					<< endl;
			mystr << "Nb streams : " << carma_host_obj_handler->get_nbStreams();
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			long *dims = carma_host_obj_handler->getDims();
			mystr << "Dims : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_COMPLEX) {
			carma_host_obj<cuDoubleComplex> *carma_host_obj_handler =
					(carma_host_obj<cuDoubleComplex> *) (handler->carma_host_object);
			mystr << "Carma Host Object : " << endl;
			mystr << "Type : double complex";
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			mystr << "Allocation : " << carma_host_obj_handler->getMetAlloc()
					<< endl;
			mystr << "Nb streams : " << carma_host_obj_handler->get_nbStreams();
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			long *dims = carma_host_obj_handler->getDims();
			mystr << "Dims : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else if (handler->type == Y_SCOMPLEX) {
			carma_host_obj<cuDoubleComplex> *carma_host_obj_handler =
					(carma_host_obj<cuDoubleComplex> *) (handler->carma_host_object);
			mystr << "Carma Host Object : " << endl;
			mystr << "Type : simple complex";
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			mystr << "Allocation : " << carma_host_obj_handler->getMetAlloc()
					<< endl;
			mystr << "Nb streams : " << carma_host_obj_handler->get_nbStreams();
			y_print(mystr.str().c_str(), 1);
			mystr.str("");
			long *dims = carma_host_obj_handler->getDims();
			mystr << "Dims : " << dims[1];
			for (int i = 2; i <= dims[0]; i++)
				mystr << "x" << dims[i];
		} else
			throw "Type unknown";
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	}
	y_print(mystr.str().c_str(), 1);
}

void yHostObj_eval(void *obj, int n)
/** @brief yHostObj_struct evaluator.
 *  @param[in] obj : yHostObj_struct to evaluate
 */
{
	yHostObj_struct *handler = (yHostObj_struct *) obj;
	try {
		if (handler->type == Y_FLOAT) {
			carma_host_obj<float> *carma_host_obj_handler =
					(carma_host_obj<float> *) (handler->carma_host_object);
			float *data = ypush_f(carma_host_obj_handler->getDims());
			carma_host_obj_handler->fill_into(data);
		} else if (handler->type == Y_DOUBLE) {
			carma_host_obj<double> *carma_host_obj_handler =
					(carma_host_obj<double> *) (handler->carma_host_object);
			double *data = ypush_d(carma_host_obj_handler->getDims());
			carma_host_obj_handler->fill_into(data);
		} else if (handler->type == Y_INT) {
			carma_host_obj<int> *carma_host_obj_handler =
					(carma_host_obj<int> *) (handler->carma_host_object);
			int *data = ypush_i(carma_host_obj_handler->getDims());
			carma_host_obj_handler->fill_into(data);
		} else if (handler->type == Y_SHORT) {
			carma_host_obj<unsigned int> *carma_host_obj_handler = (carma_host_obj<
					unsigned int> *) (handler->carma_host_object);
			unsigned int *data = (unsigned int *) ypush_i(
					carma_host_obj_handler->getDims());
			carma_host_obj_handler->fill_into(data);
		} else if (handler->type == Y_COMPLEX) {
			carma_host_obj<cuDoubleComplex> *carma_host_obj_handler =
					(carma_host_obj<cuDoubleComplex> *) (handler->carma_host_object);
			cuDoubleComplex *data = (cuDoubleComplex *) ypush_z(
					carma_host_obj_handler->getDims());
			carma_host_obj_handler->fill_into(data);
		} else if (handler->type == Y_SCOMPLEX) {
			carma_host_obj<cuComplex> *carma_host_obj_handler = (carma_host_obj<
					cuComplex> *) (handler->carma_host_object);
			long * ndims_obj = carma_host_obj_handler->getDims();
			long * ndims_data = new long[ndims_obj[0] + 2];
			ndims_data[0] = ndims_obj[0] + 1;
			ndims_data[1] = 2;
			memcpy(&ndims_data[2], &ndims_obj[1], sizeof(long) * ndims_obj[0]);
			float *data = (float *) ypush_f(ndims_data);
			carma_host_obj_handler->fill_into((float2*) data);
			delete ndims_data;
		} else
			throw "Type unknown";
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	}
}

void Y_yoga_host_obj(int argc)
/** @brief yHostObj_struct constructor.
 *  @param[in] argc : command line argument(s)
 *  as for a yoga_obj several cases are handled :
 *   - a pointer to a yoga_obj is passed as an argument : copy the pointee in a new yoga_obj
 *   - a yoga_obj is passed as an argument : create a copy
 *   - a string (type of object) and an array of dimensions (Yorick convention) is passed
 *   - an existing array is passed : create a new object an fill it with array content
 *  supported types : float, double, scomplex, complex, uint, int, long
 *  the created object is pushed on the stack
 *  several memory types are supported (using appropriate keyword) :
 *   - pagelock : paged-lock memory
 *   - zerocpy  : zero-copy memory
 *   - portable : portable memory
 *   - wricomb  : write combined memory
 *   - genepin  : generic pinned memory
 *  @see Y_yoga_obj
 */
{
	long ntot;
	long dims[Y_DIMSIZE];
	void *data_input = 0L;

	// Keywords management
	static char const *knames[6] = { "pagelock", "zerocpy", "portable",
			"wricomb", "genepin", 0 };
	static long kglobs[7];
	int kiargs[6];
	int piargs[4] = { -1, -1, -1, -1 };
	yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);

	int iarg = argc - 1;
	int parg = 0;
	while (iarg >= 0) {
		iarg = yarg_kw(iarg, kglobs, kiargs);
		if (iarg >= 0) {
			if (parg < 4)
				piargs[parg++] = iarg--;
			else
				y_error("carma_host_obj takes at most 4 positional arguments");
		}
	};
	MemAlloc mallocType = MA_MALLOC;
	if (kiargs[0] >= 0) {
#ifdef DEBUG
		cout << "allocating Page-Locked Host Memory" << endl;
#endif
		mallocType = MA_PAGELOCK;
	}
	if (kiargs[1] >= 0) {
#ifdef DEBUG
		cout << "allocating Zero-Copy Host Memory" << endl;
#endif
		mallocType = MA_ZEROCPY;
	}
	if (kiargs[2] >= 0) {
#ifdef DEBUG
		cout << "allocating Portable Host Memory" << endl;
#endif
		mallocType = MA_PORTABLE;
	}
	if (kiargs[3] >= 0) {
#ifdef DEBUG
		cout << "allocating writeCombined Host Memory" << endl;
#endif
		mallocType = MA_WRICOMB;
	}
	if (kiargs[4] >= 0) {
#ifdef DEBUG
		cout << "allocating Generic pinned Host Memory" << endl;
#endif
		mallocType = MA_GENEPIN;
	}

	try {
		int oType;
		int yType = yarg_typeid(piargs[0]);
		if (yType == Y_OPAQUE) { // Copy constructor
			yHostObj_struct *handle_obj = (yHostObj_struct *) yget_obj(
					piargs[0], &yHostObj);
			yHostObj_struct *handle = (yHostObj_struct *) ypush_obj(
					&yHostObj, sizeof(yHostObj_struct));
			handle->type = handle_obj->type;
			if (handle->type == Y_FLOAT) {
				carma_host_obj<float> *carma_host_obj_handler = (carma_host_obj<
						float> *) (handle_obj->carma_host_object);
				handle->carma_host_object = new carma_host_obj<float>(
						carma_host_obj_handler, mallocType);
			} else if (handle->type == Y_DOUBLE) {
				carma_host_obj<double> *carma_host_obj_handler = (carma_host_obj<
						double> *) (handle_obj->carma_host_object);
				handle->carma_host_object = new carma_host_obj<double>(
						carma_host_obj_handler, mallocType);
			} else if (handle->type == Y_INT) {
				carma_host_obj<int> *carma_host_obj_handler =
						(carma_host_obj<int> *) (handle_obj->carma_host_object);
				handle->carma_host_object = new carma_host_obj<int>(
						carma_host_obj_handler, mallocType);
			} else if (handle->type == Y_COMPLEX) {
				carma_host_obj<cuDoubleComplex> *carma_host_obj_handler =
						(carma_host_obj<cuDoubleComplex> *) (handle_obj->carma_host_object);
				handle->carma_host_object = new carma_host_obj<cuDoubleComplex>(
						carma_host_obj_handler, mallocType);
			} else if (handle->type == Y_SCOMPLEX) {
				carma_host_obj<cuFloatComplex> *carma_host_obj_handler =
						(carma_host_obj<cuFloatComplex> *) (handle_obj->carma_host_object);
				handle->carma_host_object = new carma_host_obj<cuFloatComplex>(
						carma_host_obj_handler, mallocType);
			}
		} else { // Standard constructor

			if (yType == Y_STRING) { // based on the input type and dimensions array
				data_input = ygeta_l(piargs[1], &ntot, dims);
				char *type_data = ygets_q(piargs[0]);

				if (strcmp(type_data, "float") == 0) {
					oType = Y_FLOAT;
				} else if (strcmp(type_data, "double") == 0) {
					oType = Y_DOUBLE;
				} else if (strcmp(type_data, "int") == 0) {
					oType = Y_INT;
				} else if (strcmp(type_data, "complex") == 0) {
					oType = Y_COMPLEX;
				} else if (strcmp(type_data, "scomplex") == 0) {
					oType = Y_SCOMPLEX;
				} else if (strcmp(type_data, "uint") == 0) {
					oType = Y_SHORT;
				}
			} else { // based on the input array
				data_input = ygeta_any(piargs[0], &ntot, dims, &oType);
			}

			yHostObj_struct *handle = (yHostObj_struct *) ypush_obj(
					&yHostObj, sizeof(yHostObj_struct));
			handle->type = oType;

			if (oType == Y_FLOAT) {
				if (yType == Y_STRING) { // based on the input type and dimensions array
					handle->carma_host_object = new carma_host_obj<float>(
							(long*) data_input, mallocType);
				} else { // based on the input array
					handle->carma_host_object = new carma_host_obj<float>(dims,
							(float*) data_input, mallocType);
				}
			} else if (oType == Y_DOUBLE) {
				if (yType == Y_STRING) {
					handle->carma_host_object = new carma_host_obj<double>(
							(long*) data_input, mallocType);
				} else {
					handle->carma_host_object = new carma_host_obj<double>(dims,
							(double*) data_input, mallocType);
				}
			} else if (oType == Y_INT) {
				if (yType == Y_STRING) {
					handle->carma_host_object = new carma_host_obj<int>(
							(long*) data_input, mallocType);
				} else {
					handle->carma_host_object = new carma_host_obj<int>(dims,
							(int*) data_input, mallocType);
				}
			} else if (oType == Y_SHORT) {
				if (yType == Y_STRING) {
					handle->carma_host_object = new carma_host_obj<unsigned int>(
							(long*) data_input, mallocType);
				} else {
					handle->carma_host_object = new carma_host_obj<unsigned int>(
							dims, (uint*) data_input, mallocType);
				}
			} else if (oType == Y_COMPLEX) {
				if (yType == Y_STRING) {
					handle->carma_host_object =
							new carma_host_obj<cuDoubleComplex>(
									(long*) data_input, mallocType);
				} else {
					handle->carma_host_object =
							new carma_host_obj<cuDoubleComplex>(dims,
									(cuDoubleComplex*) data_input, mallocType);
				}
			} else if (oType == Y_SCOMPLEX) {
				if (yType == Y_STRING) {
					handle->carma_host_object =
							new carma_host_obj<cuFloatComplex>(
									(long*) data_input, mallocType);
				} else {
					stringstream buf;
					buf << "Type not supported in yorick in " << __FILE__ << "@"
							<< __LINE__ << endl;
					throw buf.str();
				}
			} else {
				stringstream buf;
				buf << "Type not found in " << __FILE__ << "@" << __LINE__
						<< endl;
				throw buf.str();
			}
		}
	} catch (string &msg) {
		y_error(msg.c_str());
	} catch (char const * msg) {
		y_error(msg);
	} catch (...) {
		stringstream buf;
		buf << "unknown error with yObj construction in " << __FILE__ << "@"
				<< __LINE__ << endl;
		y_error(buf.str().c_str());
	}
}

void Y_yoga_host_getp(int argc)
/** @brief simple routine to retrieve the address to a yoga_host_obj
 *  @param[in] argc : command line argument (yoga_host_obj expected)
 *  the corresponding address is pushed on the stack
 */
{
	yHostObj_struct *handle = (yHostObj_struct *) yget_obj(argc - 1,
			&yHostObj);
	long a = 0;
	ypointer_t *ptr_handle = ypush_p(&a);
	int type = handle->type;
	if (type == Y_DOUBLE) {
		carma_host_obj<double>* test =
				(carma_host_obj<double>*) handle->carma_host_object;
		*ptr_handle = test->getData();
	}
}

//#ifdef _USE_MAGMA
/*
 *  _ __ ___   __ _  __ _ _ __ ___   __
 * | '_ ` _ \ / _` |/ _` | '_ ` _ \ / _` |
 * | | | | | | (_| | (_| | | | | | | (_| |
 * |_| |_| |_|\__,_|\__, |_| |_| |_|\__,_|
 *                  |___/
 */

void Y_yoga_svd(int argc)
/** @brief wrapper routine for yoga svd method using magma
 *  @param[in] argc : command line arguments
 *    - first  : the yoga_obj to be decomposed
 *    - second : a yoga_obj for eignevalues
 *    - third  : a yoga_obj for thr U matrix
 *    - fourth : a yoga_obj for the VT matrix
 *  only floating point matrices can be decomposed (single or double precision)
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle_mat = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_eigenvals = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		yObj_struct *handle_U = (yObj_struct *) yget_obj(argc - 3, &yObj);
		yObj_struct *handle_VT = (yObj_struct *) yget_obj(argc - 4, &yObj);

		if (handle_mat->type == Y_FLOAT) {
			caObjS *carma_obj_handler_mat = (caObjS *) (handle_mat->carma_object);
			caObjS *carma_obj_handler_eigenvals =
					(caObjS *) (handle_eigenvals->carma_object);
			caObjS *carma_obj_handler_U = (caObjS *) (handle_U->carma_object);
			caObjS *carma_obj_handler_VT = (caObjS *) (handle_VT->carma_object);
			carma_svd(carma_obj_handler_mat, carma_obj_handler_eigenvals,
					carma_obj_handler_VT, carma_obj_handler_U);
		}
		if (handle_mat->type == Y_DOUBLE) {
			caObjD *carma_obj_handler_mat = (caObjD *) (handle_mat->carma_object);
			caObjD *carma_obj_handler_eigenvals =
					(caObjD *) (handle_eigenvals->carma_object);
			caObjD *carma_obj_handler_U = (caObjD *) (handle_U->carma_object);
			caObjD *carma_obj_handler_VT = (caObjD *) (handle_VT->carma_object);
			carma_svd(carma_obj_handler_mat, carma_obj_handler_eigenvals,
					carma_obj_handler_VT, carma_obj_handler_U);
		}
	} else {
	}
}

void Y_yoga_svd_host(int argc)
/** @brief wrapper routine for yoga svd method using magma on a yoga_host_obj
 *  @param[in] argc : command line arguments
 *    - first  : the yoga_host_obj to be decomposed
 *    - second : a yoga_host_obj for eignevalues
 *    - third  : a yoga_host_obj for thr U matrix
 *    - fourth : a yoga_host_obj for the VT matrix
 *  only floating point matrices can be decomposed (single or double precision)
 */
{
	if (yarg_subroutine()) {
		yHostObj_struct *handle_mat = (yHostObj_struct *) yget_obj(argc - 1,
				&yHostObj);
		yHostObj_struct *handle_eigenvals = (yHostObj_struct *) yget_obj(
				argc - 2, &yHostObj);
		yHostObj_struct *handle_U = (yHostObj_struct *) yget_obj(argc - 3,
				&yHostObj);
		yHostObj_struct *handle_VT = (yHostObj_struct *) yget_obj(argc - 4,
				&yHostObj);

		if (handle_mat->type == Y_FLOAT) {
			carma_host_obj<float> *carma_obj_handler_mat =
					(carma_host_obj<float> *) (handle_mat->carma_host_object);

			carma_host_obj<float> *carma_obj_handler_eigenvals = (carma_host_obj<
					float> *) (handle_eigenvals->carma_host_object);

			carma_host_obj<float> *carma_obj_handler_U =
					(carma_host_obj<float> *) (handle_U->carma_host_object);

			carma_host_obj<float> *carma_obj_handler_VT =
					(carma_host_obj<float> *) (handle_VT->carma_host_object);

			carma_svd(carma_obj_handler_mat, carma_obj_handler_eigenvals,
					carma_obj_handler_VT, carma_obj_handler_U);
		}
		if (handle_mat->type == Y_DOUBLE) {
			carma_host_obj<double> *carma_obj_handler_mat =
					(carma_host_obj<double> *) (handle_mat->carma_host_object);

			carma_host_obj<double> *carma_obj_handler_eigenvals = (carma_host_obj<
					double> *) (handle_eigenvals->carma_host_object);

			carma_host_obj<double> *carma_obj_handler_U =
					(carma_host_obj<double> *) (handle_U->carma_host_object);

			carma_host_obj<double> *carma_obj_handler_VT =
					(carma_host_obj<double> *) (handle_VT->carma_host_object);

			carma_svd(carma_obj_handler_mat, carma_obj_handler_eigenvals,
					carma_obj_handler_VT, carma_obj_handler_U);
		}
	} else {
	}
}
//#endif

/*
 *             _
 *   ___ _   _| | __ _
 *  / __| | | | |/ _` |
 * | (__| |_| | | (_| |
 *  \___|\__,_|_|\__,_|
 */

void Y_yoga_cula_svd(int argc)
/** @brief wrapper routine for yoga svd method using cula
 *  @param[in] argc : command line arguments
 *    - first  : the yoga_obj to be decomposed
 *    - second : a yoga_obj for eignevalues
 *    - third  : a yoga_obj for thr U matrix
 *    - fourth : a yoga_obj for the VT matrix
 *  only floating point matrices can be decomposed (single or double precision)
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle_mat = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_eigenvals = (yObj_struct *) yget_obj(argc - 2,
				&yObj);
		yObj_struct *handle_U = (yObj_struct *) yget_obj(argc - 3, &yObj);
		yObj_struct *handle_VT = (yObj_struct *) yget_obj(argc - 4, &yObj);

		if (handle_mat->type == Y_FLOAT) {
			caObjS *carma_obj_handler_mat = (caObjS *) (handle_mat->carma_object);
			caObjS *carma_obj_handler_eigenvals =
					(caObjS *) (handle_eigenvals->carma_object);
			caObjS *carma_obj_handler_U = (caObjS *) (handle_U->carma_object);
			caObjS *carma_obj_handler_VT = (caObjS *) (handle_VT->carma_object);
			carma_cula_svd(carma_obj_handler_mat, carma_obj_handler_eigenvals,
					carma_obj_handler_VT, carma_obj_handler_U);
		}
		if (handle_mat->type == Y_DOUBLE) {
			caObjD *carma_obj_handler_mat = (caObjD *) (handle_mat->carma_object);
			caObjD *carma_obj_handler_eigenvals =
					(caObjD *) (handle_eigenvals->carma_object);
			caObjD *carma_obj_handler_U = (caObjD *) (handle_U->carma_object);
			caObjD *carma_obj_handler_VT = (caObjD *) (handle_VT->carma_object);
			carma_cula_svd(carma_obj_handler_mat, carma_obj_handler_eigenvals,
					carma_obj_handler_VT, carma_obj_handler_U);
		}
	} else {
	}
}

void Y_yoga_cula_svd_host(int argc)
/** @brief wrapper routine for yoga svd method using cula on a yoga_host_obj
 *  @param[in] argc : command line arguments
 *    - first  : the yoga_host_obj to be decomposed
 *    - second : a yoga_host_obj for eignevalues
 *    - third  : a yoga_host_obj for thr U matrix
 *    - fourth : a yoga_host_obj for the VT matrix
 *  only floating point matrices can be decomposed (single or double precision)
 */
{
	if (yarg_subroutine()) {
		yHostObj_struct *handle_mat = (yHostObj_struct *) yget_obj(argc - 1,
				&yHostObj);
		yHostObj_struct *handle_eigenvals = (yHostObj_struct *) yget_obj(
				argc - 2, &yHostObj);
		yHostObj_struct *handle_U = (yHostObj_struct *) yget_obj(argc - 3,
				&yHostObj);
		yHostObj_struct *handle_VT = (yHostObj_struct *) yget_obj(argc - 4,
				&yHostObj);

		if (handle_mat->type == Y_FLOAT) {
			carma_host_obj<float> *carma_obj_handler_mat =
					(carma_host_obj<float> *) (handle_mat->carma_host_object);

			carma_host_obj<float> *carma_obj_handler_eigenvals = (carma_host_obj<
					float> *) (handle_eigenvals->carma_host_object);

			carma_host_obj<float> *carma_obj_handler_U =
					(carma_host_obj<float> *) (handle_U->carma_host_object);

			carma_host_obj<float> *carma_obj_handler_VT =
					(carma_host_obj<float> *) (handle_VT->carma_host_object);

			carma_cula_svd(carma_obj_handler_mat, carma_obj_handler_eigenvals,
					carma_obj_handler_VT, carma_obj_handler_U);
		}
		if (handle_mat->type == Y_DOUBLE) {
			carma_host_obj<double> *carma_obj_handler_mat =
					(carma_host_obj<double> *) (handle_mat->carma_host_object);

			carma_host_obj<double> *carma_obj_handler_eigenvals = (carma_host_obj<
					double> *) (handle_eigenvals->carma_host_object);

			carma_host_obj<double> *carma_obj_handler_U =
					(carma_host_obj<double> *) (handle_U->carma_host_object);

			carma_host_obj<double> *carma_obj_handler_VT =
					(carma_host_obj<double> *) (handle_VT->carma_host_object);

			carma_cula_svd(carma_obj_handler_mat, carma_obj_handler_eigenvals,
					carma_obj_handler_VT, carma_obj_handler_U);
		}
	} else {
	}
}

/*
 *
 *    __ _ _ __ _ __ __ _ _   _ ___
 *   / _` | '__| '__/ _` | | | / __|
 *  | (_| | |  | | | (_| | |_| \__ \
 *   \__,_|_|  |_|  \__,_|\__, |___/
 *                         |___/
 */

void Y_yoga_getarray(int argc)
/*
 subroutine or function accepted
 interprets void as the full range of the corresponding dimension
 */
{
	yObj_struct *handle_out = NULL;
	yObj_struct *handle_in = NULL;

	if (yarg_subroutine()) {
		handle_out = (yObj_struct *) yget_obj(argc - 1, &yObj);
		handle_in = (yObj_struct *) yget_obj(argc - 2, &yObj);
		if (handle_in->device != handle_out->device)
			y_error("getarray only on the same device");
		argc--;
	} else {
		handle_in = (yObj_struct *) yget_obj(argc - 1, &yObj);
	}
	carma_context *context_handle = _getCurrentContext();
	context_handle->set_activeDeviceForCpy(handle_in->device);

	if ((yarg_typeid(argc - 2) != Y_RANGE) && (yarg_typeid(argc - 2) != Y_VOID))
		y_error("expecting a range");
	long mrange[8];
	int Nlig;
	int Ncol;
	if (yarg_typeid(argc - 2) != Y_VOID)
		mrange[0] = yget_range(argc - 2, &mrange[1]);
	else
		mrange[0] = mrange[1] = mrange[2] = mrange[3] = 0;
	if (argc > 2) {
		if ((yarg_typeid(argc - 3) != Y_RANGE)
				&& (yarg_typeid(argc - 3) != Y_VOID))
			y_error("expecting a range");
		if (yarg_typeid(argc - 3) != Y_VOID)
			mrange[4] = yget_range(argc - 3, &mrange[5]);
		else
			mrange[4] = mrange[5] = mrange[6] = mrange[7] = 0;
	} else
		mrange[4] = mrange[5] = mrange[6] = mrange[7] = 0;

	caObjS *carma_in_handler = (caObjS *) (handle_in->carma_object);
	caObjS *carma_out_handler = NULL;

	if (mrange[1] == 0) {
		mrange[1] = 1;
		mrange[2] = carma_in_handler->getDims(1);
	} else if ((mrange[2] - mrange[1] + 1 > carma_in_handler->getDims(1))
			|| mrange[2] > carma_in_handler->getDims(1))
		y_error("range out of bounds");
	if (mrange[4] == 0) {
		mrange[5] = 1;
		mrange[6] = carma_in_handler->getDims(2);
	} else if ((mrange[6] - mrange[5] + 1 > carma_in_handler->getDims(2))
			|| mrange[6] > carma_in_handler->getDims(2))
		y_error("range out of bounds");
	Ncol = mrange[2] - mrange[1] + 1;
	Nlig = mrange[6] - mrange[5] + 1;
	long dims[3];
	dims[0] = 2;
	dims[1] = Ncol;
	dims[2] = Nlig;

	if (!yarg_subroutine()) {
		handle_out = (yObj_struct *) ypush_obj(&yObj, sizeof(yObj_struct));
		handle_out->device = handle_in->device;
		if (handle_in->type == Y_FLOAT) {
			handle_out->type = Y_FLOAT;
			handle_out->carma_object = new caObjS(context_handle, dims);
			carma_out_handler = (caObjS *) (handle_out->carma_object);
		}
	} else {
		carma_out_handler = (caObjS *) (handle_out->carma_object);
	}

	int x0 = carma_in_handler->getDims(1) * (mrange[5] - 1) + (mrange[2] - 1);
	//if (Ncol == 1) x0 = (mrange[1]-1);
	getarray2d(carma_out_handler->getData(), carma_in_handler->getData(), x0,
			Ncol, carma_in_handler->getDims(1), Nlig * Ncol);
}

void Y_yoga_fillarray(int argc)
/*
 fills the content of the first arg with the content of the second args
 in the range specified by third and fourth args
 only subroutine accepted
 interprets void as the full range of the corresponding dimension
 */
{
	if (yarg_subroutine()) {
		yObj_struct *handle_out = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_in = (yObj_struct *) yget_obj(argc - 2, &yObj);
		if (handle_in->device != handle_out->device)
			y_error("getarray only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDeviceForCpy(handle_in->device);
		caObjS *carma_out_handler = (caObjS *) (handle_out->carma_object);
		caObjS *carma_in_handler = (caObjS *) (handle_in->carma_object);
		argc--;
		if ((yarg_typeid(argc - 2) != Y_RANGE)
				&& (yarg_typeid(argc - 2) != Y_VOID))
			y_error("expecting a range");
		long mrange[8];
		int Nlig;
		int Ncol;
		if (yarg_typeid(argc - 2) != Y_VOID)
			mrange[0] = yget_range(argc - 2, &mrange[1]);
		else
			mrange[0] = mrange[1] = mrange[2] = mrange[3] = 0;
		if (argc > 2) {
			if ((yarg_typeid(argc - 3) != Y_RANGE)
					&& (yarg_typeid(argc - 3) != Y_VOID))
				y_error("expecting a range");
			if (yarg_typeid(argc - 3) != Y_VOID)
				mrange[4] = yget_range(argc - 3, &mrange[5]);
			else
				mrange[4] = mrange[5] = mrange[6] = mrange[7] = 0;
		} else
			mrange[4] = mrange[5] = mrange[6] = mrange[7] = 0;

		if (mrange[1] == 0) {
			mrange[1] = 1;
			mrange[2] = carma_out_handler->getDims(1);
		} else if ((mrange[2] - mrange[1] + 1 > carma_out_handler->getDims(1))
				|| mrange[2] > carma_out_handler->getDims(1))
			y_error("range out of bounds");
		if (mrange[4] == 0) {
			mrange[5] = 1;
			mrange[6] = carma_out_handler->getDims(2);
		} else if ((mrange[6] - mrange[5] + 1 > carma_out_handler->getDims(2))
				|| mrange[6] > carma_out_handler->getDims(2))
			y_error("range out of bounds");
		Ncol = mrange[2] - mrange[1] + 1;
		Nlig = mrange[6] - mrange[5] + 1;
		long dims[3];
		dims[0] = 2;
		dims[1] = Ncol;
		dims[2] = Nlig;

		int x0 = carma_out_handler->getDims(1) * (mrange[5] - 1)
				+ (mrange[1] - 1);
		//if (Ncol == 1) x0 = (mrange[1]-1);
		fillarray2d(carma_out_handler->getData(), carma_in_handler->getData(), x0,
				Ncol, carma_out_handler->getDims(1), Nlig * Ncol);
	} else
		y_error("can only be called as a subroutine");
}

void Y_yoga_getvalue(int argc)
/*
 useful to test array indexing
 on the gpu the linear indexing is column / line meaning :
 element #0 is the 1rst elem of the array (1rst of 1rst line and 1rst of 1rst col
 element #1 is the 2nd elem of 1rst line and 1rst elem of 2nd column
 element #2 is the 3rd elem of 1rst line and 1rst elem of 3rd colum
 elem # ncol is the 1rst elem of 2nd line and 2nd elem of 1rst colum
 etc ...
 in yorick the same applies
 in 2d (yorick only)
 elem (1,2) is elem # ncol
 elem (2,1) is elem #
 so the translation is done as : # line : (Nelem / ncol)
 # col  : (Nelem % ncol)
 */
{
	if (yarg_subroutine())
		y_error("can only be called as a function");
	else {
		yObj_struct *handle_in = (yObj_struct *) yget_obj(argc - 1, &yObj);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDeviceForCpy(handle_in->device);
		int pos = ygets_i(argc-2);
		long dims[2];
		dims[0] = dims[1] = 1;
		caObjS *carma_in_handler = (caObjS *) (handle_in->carma_object);
		float res;
		cudaMemcpy(&res, &((carma_in_handler->getData())[pos - 1]),
				sizeof(float), cudaMemcpyDeviceToHost);
		float* yres = ypush_f(dims);
		yres[0] = res;
	}
}

void Y_yoga_plus(int argc) {
	if (yarg_subroutine()) {
		yObj_struct *handle_in = (yObj_struct *) yget_obj(argc - 1, &yObj);
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_in->device);
		if (handle_in->type == Y_FLOAT) {
			float alpha = ygets_f(argc-2);
			caObjS *carma_in_handler = (caObjS *) (handle_in->carma_object);
			carma_plus(carma_in_handler->getData(), alpha,
					carma_in_handler->getNbElem());
		} else
			y_error("double not implemented yet");
	} else
		y_error("can only be called as a function");
}

void Y_yoga_plusai(int argc) {
	if (yarg_subroutine()) {
		yObj_struct *handle_out = (yObj_struct *) yget_obj(argc - 1,
				&yObj);
		yObj_struct *handle_in = (yObj_struct *) yget_obj(argc - 2, &yObj);
		if (handle_in->device != handle_out->device)
			y_error("getarray only on the same device");
		carma_context *context_handle = _getCurrentContext();
		context_handle->set_activeDevice(handle_in->device);
		int idx = ygets_i(argc-3);
		int sgn;
		if (argc > 3)
			sgn = ygets_i(argc-4);
		else
			sgn = 1;
		if ((sgn != -1) && (sgn != 1))
			y_error("arg 4 must be +1 or -1");
		if (handle_out->type == Y_FLOAT) {
			caObjS *carma_in_handler = (caObjS *) (handle_in->carma_object);
			caObjS *carma_out_handler = (caObjS *) (handle_out->carma_object);
			carma_plusai(carma_out_handler->getData(), carma_in_handler->getData(),
					idx, sgn, carma_out_handler->getNbElem());
		} else
			y_error("double not implemented yet");
	} else
		y_error("can only be called as a function");
}

void Y_yoga_test(int argc) {
	if (yarg_subroutine()) {
		int test = ygets_i(argc-1);
		fprintf(stderr, "%d\n", test);
	}
}

yObj_struct* yoga_getyObj(int argc, int pos) {
	return (yObj_struct *) yget_obj(argc - pos, &yObj);
}

yObj_struct* yoga_pushyObj(void) {
	return (yObj_struct *) ypush_obj(&yObj, sizeof(yObj_struct));
}

}
