#include "convolutionFFT2D_common.h"
#include <string>
#include <carma_obj.h>

////////////////////////////////////////////////////////////////////////////////
// Helper functions
/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */
////////////////////////////////////////////////////////////////////////////////
int snapTransformSize(unsigned int dataSize) {
	int hiBit;
	unsigned int lowPOT, hiPOT;

	dataSize = iAlignUp(dataSize, 16);

	for (hiBit = 31; hiBit >= 0; hiBit--)
		if (dataSize & (1U << hiBit))
			break;

	lowPOT = 1U << hiBit;
	if (lowPOT == dataSize)
		return dataSize;

	hiPOT = 1U << (hiBit + 1);
	if (hiPOT <= 1024)
		return hiPOT;
	else
		return iAlignUp(dataSize, 512);
}
/*
 ____                      _           _____ _____ _____ 
 / ___|___  _ ____   _____ | |_   _____|  ___|  ___|_   _|
 | |   / _ \| '_ \ \ / / _ \| \ \ / / _ \ |_  | |_    | |  
 | |__| (_) | | | \ V / (_) | |\ V /  __/  _| |  _|   | |  
 \____\___/|_| |_|\_/ \___/|_| \_/ \___|_|   |_|     |_|  
 
 */

int carma_initfftconv(caObjS *data_in, caObjS *kernel_in, caObjS *padded_data,
		caObjC *padded_spectrum, int kernelY, int kernelX) {

	float *odata = padded_data->getOData();
	cufftHandle *plan = padded_data->getPlan(); ///< FFT plan

	if (odata == NULL)
		cutilSafeCall(
				cudaMalloc((void** )&odata,
						sizeof(float) * padded_data->getNbElem()));
	cutilSafeCall(
			cudaMemset(padded_data->getOData(), 0,
					padded_data->getNbElem() * sizeof(float)));
	if (data_in->getDims(0) == 3) {
		int mdims[2];
		mdims[0] = padded_data->getDims(1);
		mdims[1] = padded_data->getDims(2);
		cufftSafeCall(
				cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_R2C ,padded_data->getDims(3)));
	} else {
		cufftSafeCall(
				cufftPlan2d(plan, padded_data->getDims(1),
						padded_data->getDims(2), CUFFT_R2C));
	}
	// plan fwd
	cuFloatComplex *odata2 = padded_spectrum->getOData();
	if (odata2 == NULL)
		cutilSafeCall(
				cudaMalloc((void** )&odata2,
						sizeof(cuFloatComplex) * padded_spectrum->getNbElem()));
	cutilSafeCall(
			cudaMemset(padded_spectrum->getOData(), 0,
					padded_spectrum->getNbElem() * sizeof(float)));

	plan = padded_spectrum->getPlan();
	if (data_in->getDims(0) == 3) {
		int mdims[2];
		mdims[0] = (int) padded_data->getDims(1);
		mdims[1] = (int) padded_data->getDims(2);
		cufftSafeCall(
				cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2R ,(int)padded_spectrum->getDims(3)));
	} else {
		cufftSafeCall(
				cufftPlan2d(plan, padded_data->getDims(1),
						padded_data->getDims(2), CUFFT_C2R));
	}
	// plan inv

	if (kernel_in->getDims(0) == 3) {
		padKernel3d(padded_data->getOData(), kernel_in->getData(),
				padded_data->getDims(1), padded_data->getDims(2),
				kernel_in->getDims(1), kernel_in->getDims(2), kernelY, kernelX,
				kernel_in->getDims(3));
	} else {
		padKernel(padded_data->getOData(), kernel_in->getData(),
				padded_data->getDims(1), padded_data->getDims(2),
				kernel_in->getDims(1), kernel_in->getDims(2), kernelY, kernelX);
	}

	if (data_in->getDims(0) == 3) {
		padDataClampToBorder3d(padded_data->getData(), data_in->getData(),
				padded_data->getDims(1), padded_data->getDims(2),
				data_in->getDims(1), data_in->getDims(2), kernel_in->getDims(1),
				kernel_in->getDims(2), kernelY, kernelX, data_in->getDims(3));
	} else {
		padDataClampToBorder(padded_data->getData(), data_in->getData(),
				padded_data->getDims(1), padded_data->getDims(2),
				data_in->getDims(1), data_in->getDims(2), kernel_in->getDims(1),
				kernel_in->getDims(2), kernelY, kernelX);
	}

	return EXIT_SUCCESS;
}

int carma_fftconv(caObjS *data_out, caObjS *padded_data,
		caObjC *padded_spectrum, int kernelY, int kernelX) {

	carma_fft(padded_data->getOData(), padded_spectrum->getOData(), 1,
			*padded_data->getPlan());
	// keyword dir in fft is not relevant in this case
	carma_fft(padded_data->getData(), padded_spectrum->getData(), 1,
			*padded_data->getPlan());

	int nim;
	nim = data_out->getDims(0) == 3 ? data_out->getDims(3) : 1;

	modulateAndNormalize((fComplex *) (padded_spectrum->getData()),
			(fComplex *) (padded_spectrum->getOData()), padded_data->getDims(1),
			padded_data->getDims(2), 1, nim);

	carma_fft(padded_spectrum->getData(), padded_data->getData(), 1,
			*padded_spectrum->getPlan());

	int N = padded_data->getDims(2) * padded_data->getDims(1);
	int n = data_out->getDims(1) * data_out->getDims(2);

	fftconv_unpad(data_out->getData(), padded_data->getData(),
			padded_data->getDims(2), data_out->getDims(1), data_out->getDims(2),
			N, n, nim);

	return EXIT_SUCCESS;
}

