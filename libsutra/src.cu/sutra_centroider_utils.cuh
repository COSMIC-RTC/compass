#ifndef SUTRA_CENTROIDER_UTILS_CUH
#define SUTRA_CENTROIDER_UTILS_CUH

#include <sutra_centroider.hpp>

namespace sutra
{

    struct SlopesIndex {
        int32_t factor, offset;

        __host__ __device__ SlopesIndex(int32_t nvalid, SlopeOrder so):
            factor((so == SlopeOrder::untied) ? 1 : 2),
            offset((so == SlopeOrder::untied) ? nvalid : 1)
        {}

        __host__ __device__
        constexpr int32_t x(std::size_t pos) const noexcept
        {
            return factor * pos;
        }

        __host__ __device__
        constexpr int32_t y(std::size_t pos) const noexcept
        {
            return factor * pos + offset;
        }

    };


} // namespace sutra

#endif //SUTRA_CENTROIDER_UTILS_CUH