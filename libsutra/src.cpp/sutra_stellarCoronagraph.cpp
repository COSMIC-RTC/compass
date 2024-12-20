// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_stellarCoronagraph.cpp
//! \ingroup   libsutra
//! \class     SutraStellarCoronagraph
//! \brief     this class provides the Stellar coronograph features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_stellarCoronagraph.hpp>

SutraStellarCoronagraph::SutraStellarCoronagraph(CarmaContext *context, SutraSource *d_source,
                                    int32_t im_dimx, int32_t im_dimy, int32_t fpm_dimx, int32_t fpm_dimy,
                                    float *wavelength, int32_t nWavelength, bool babinet,
                                    int32_t device):
    SutraCoronagraph(context, "perfect", d_source, im_dimx, im_dimy, wavelength, nWavelength, device),
    fpmDimx(fpm_dimx),
    fpmDimy(fpm_dimy),
    babinet(babinet) {
        AA = {
            {"img", make_tuple(mftVec(nWavelength, nullptr), vector<int64_t>({2, imageDimx, pupDimx}))},
            {"psf", make_tuple(mftVec(nWavelength, nullptr), vector<int64_t>({2, imageDimx, pupDimx}))},
            {"fpm", make_tuple(mftVec(nWavelength, nullptr), vector<int64_t>({2, fpmDimx, pupDimx}))},
            {"lyot", make_tuple(mftVec(nWavelength, nullptr), vector<int64_t>({2, pupDimx, fpmDimx}))}
        };
        BB = {
            {"img", make_tuple(mftVec(nWavelength, nullptr), vector<int64_t>({2, pupDimy, imageDimy}))},
            {"psf", make_tuple(mftVec(nWavelength, nullptr), vector<int64_t>({2, pupDimy, imageDimy}))},
            {"fpm", make_tuple(mftVec(nWavelength, nullptr), vector<int64_t>({2, pupDimy, fpmDimy}))},
            {"lyot", make_tuple(mftVec(nWavelength, nullptr), vector<int64_t>({2, fpmDimy, pupDimy}))}
        };
        norm = {
            {"img", vector<float>(nWavelength, 1)},
            {"psf", vector<float>(nWavelength, 1)},
            {"fpm", vector<float>(nWavelength, 1)},
            {"lyot", vector<float>(nWavelength, 1)}
        };
        tmp_mft = {
            {"img", make_tuple(nullptr, vector<int64_t>({2, imageDimx, pupDimy}))},
            {"psf", make_tuple(nullptr, vector<int64_t>({2, imageDimx, pupDimy}))},
            {"fpm", make_tuple(nullptr, vector<int64_t>({2, fpmDimx, pupDimy}))},
            {"lyot", make_tuple(nullptr, vector<int64_t>({2, pupDimx, fpmDimy}))}
        };

        int64_t dims[3] = {2, pupDimx, pupDimy};
        d_apodizer = new CarmaObj<float>(current_context, dims);
        d_lyot_stop = new CarmaObj<float>(current_context, dims);
        d_electric_field_babinet = new CarmaObj<cuFloatComplex>(current_context, dims);
        dims[1] = fpmDimx;
        dims[2] = fpmDimy;
        d_electric_field_fpm = new CarmaObj<cuFloatComplex>(current_context, dims);
}

int32_t SutraStellarCoronagraph::set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm0,
                                        std::string mftType) {
    if(AA.count(mftType) < 1) {
        std::cout << "Invalid mftType. Must be img or psf" << std::endl;
        return EXIT_FAILURE;
    }

    vector<int64_t> dims;
    CarmaObj<cuFloatComplex> *data;

    for (int32_t i = 0; i < wavelength.size() ; i++) {
        dims = std::get<1>(AA[mftType]);
        std::get<0>(AA[mftType])[i] = new CarmaObj<cuFloatComplex>(current_context, dims.data(), A + i * dims[1] * dims[2]);
        dims = std::get<1>(BB[mftType]);
        std::get<0>(BB[mftType])[i] = new CarmaObj<cuFloatComplex>(current_context, dims.data(), B + i * dims[1] * dims[2]);
        norm[mftType][i] = norm0[i];
    }

    if (std::get<0>(tmp_mft[mftType]) == nullptr) {
        dims = std::get<1>(tmp_mft[mftType]);
        std::get<0>(tmp_mft[mftType]) = new CarmaObj<cuFloatComplex>(current_context, dims.data());
    }
    return EXIT_SUCCESS;
}

int32_t SutraStellarCoronagraph::set_focal_plane_mask(float *mask) {
    int64_t dims[3] = {2, fpmDimx, fpmDimy};

    if (focal_plane_mask.empty()) {
        for (int32_t i = 0; i < wavelength.size() ; i++) {
            focal_plane_mask.push_back(new CarmaObj<float>(current_context, dims,
                                                            mask + i * dims[1] * dims[2]));
        }
    }
    else {
        for (int32_t i = 0; i < wavelength.size() ; i++) {
            focal_plane_mask[i]->host2device(mask + i * dims[1] * dims[2]);
        }
    }

    return EXIT_SUCCESS;
}

int32_t SutraStellarCoronagraph::set_apodizer(float *apodizer) {
    d_apodizer->host2device(apodizer);
    return EXIT_SUCCESS;
}

int32_t SutraStellarCoronagraph::set_lyot_stop(float *lyot_stop) {
    d_lyot_stop->host2device(lyot_stop);
    return EXIT_SUCCESS;
}

int32_t SutraStellarCoronagraph::_compute_image(bool center_on_pixel, bool accumulate, bool no_fpm) {
    std::vector<CarmaObj<cuFloatComplex>*> mftA = std::get<0>(AA["img"]);
    std::vector<CarmaObj<cuFloatComplex>*> mftB = std::get<0>(BB["img"]);
    CarmaObj<cuFloatComplex>* tmp = std::get<0>(tmp_mft["img"]);
    std::vector<float> mftNorm = norm["img"];
    if (center_on_pixel) {
        mftA = std::get<0>(AA["psf"]);
        mftB = std::get<0>(BB["psf"]);
        tmp = std::get<0>(tmp_mft["psf"]);
        mftNorm = norm["psf"];
    }

    d_image_se->reset();
    for (int32_t i = 0 ; i < wavelength.size(); i++) {
        compute_electric_field(i);
        apply_mask(d_electric_field->get_data(), d_apodizer->get_data(),
                    d_electric_field->get_nb_elements(), current_context->get_device(device));
        mft(std::get<0>(AA["fpm"])[i], std::get<0>(BB["fpm"])[i], std::get<0>(tmp_mft["fpm"]),
            d_electric_field, d_electric_field_fpm, norm["fpm"][i]);

        if(!no_fpm)
            apply_mask(d_electric_field_fpm->get_data(), focal_plane_mask[i]->get_data(),
                    d_electric_field_fpm->get_nb_elements(), current_context->get_device(device));

        if (babinet && !no_fpm) {
            mft(std::get<0>(AA["lyot"])[i], std::get<0>(BB["lyot"])[i], std::get<0>(tmp_mft["lyot"]),
                d_electric_field_fpm, d_electric_field_babinet, norm["lyot"][i]);
            cuFloatComplex alpha(make_float2(-1, 0));
            d_electric_field->axpy(alpha, d_electric_field_babinet, 1, 1);
        }
        if(!babinet) {
            mft(std::get<0>(AA["lyot"])[i], std::get<0>(BB["lyot"])[i], std::get<0>(tmp_mft["lyot"]),
                d_electric_field_fpm, d_electric_field, norm["lyot"][i]);
        }
        apply_mask(d_electric_field->get_data(), d_lyot_stop->get_data(),
                    d_electric_field->get_nb_elements(), current_context->get_device(device));
        mft(mftA[i], mftB[i], tmp, d_electric_field, d_complex_image, mftNorm[i]);

        accumulate_abs2(d_complex_image->get_data(), d_image_se->get_data(),
                        d_image_se->get_nb_elements(), current_context->get_device(device));
    }

    if(accumulate) {
        d_image_le->axpy(1.0f, d_image_se, 1, 1);
        cntImg += 1;
    }
    return EXIT_SUCCESS;
}

int32_t SutraStellarCoronagraph::compute_psf(bool accumulate) {
    d_psf_se->reset();
    for (int32_t i = 0 ; i < wavelength.size(); i++) {
        compute_electric_field(i);
        apply_mask(d_electric_field->get_data(), d_apodizer->get_data(),
                    d_electric_field->get_nb_elements(), current_context->get_device(device));
        apply_mask(d_electric_field->get_data(), d_lyot_stop->get_data(),
                    d_electric_field->get_nb_elements(), current_context->get_device(device));
        mft(std::get<0>(AA["psf"])[i], std::get<0>(BB["psf"])[i], std::get<0>(tmp_mft["psf"]), d_electric_field, d_complex_image, norm["psf"][i]);

        accumulate_abs2(d_complex_image->get_data(), d_psf_se->get_data(),
                        d_psf_se->get_nb_elements(), current_context->get_device(device));
    }

    if(accumulate) {
        d_psf_le->axpy(1.0f, d_psf_se, 1, 1);
        cntPsf += 1;
    }

    return EXIT_SUCCESS;
}

int32_t SutraStellarCoronagraph::compute_image(bool accumulate) {
    return _compute_image(false, accumulate, false);
}

int32_t SutraStellarCoronagraph::compute_image_normalization() {
    return _compute_image(true, false, true);
}