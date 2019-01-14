################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables
CPP_SRCS += \
../src.cpp/carma_context.cpp \
../src.cpp/carma_cublas.cpp \
../src.cpp/carma_cula.cpp \
../src.cpp/carma_cusparse.cpp \
../src.cpp/carma_fft.cpp \
../src.cpp/carma_fft_conv.cpp \
../src.cpp/carma_host_obj.cpp \
../src.cpp/carma_ipcs.cpp \
../src.cpp/carma_ksparse.cpp \
../src.cpp/carma_magma.cpp \
../src.cpp/carma_multithread.cpp \
../src.cpp/carma_obj.cpp \
../src.cpp/carma_rng.cpp \
../src.cpp/carma_sparse_host_obj.cpp \
../src.cpp/carma_sparse_obj.cpp \
../src.cpp/carma_streams.cpp \
../src.cpp/carma_utils.cpp

OBJS += \
./src.cpp/carma_context.o \
./src.cpp/carma_cublas.o \
./src.cpp/carma_cula.o \
./src.cpp/carma_cusparse.o \
./src.cpp/carma_fft.o \
./src.cpp/carma_fft_conv.o \
./src.cpp/carma_host_obj.o \
./src.cpp/carma_ipcs.o \
./src.cpp/carma_ksparse.o \
./src.cpp/carma_magma.o \
./src.cpp/carma_multithread.o \
./src.cpp/carma_obj.o \
./src.cpp/carma_rng.o \
./src.cpp/carma_sparse_host_obj.o \
./src.cpp/carma_sparse_obj.o \
./src.cpp/carma_streams.o \
./src.cpp/carma_utils.o

CPP_DEPS += \
./src.cpp/carma_context.d \
./src.cpp/carma_cublas.d \
./src.cpp/carma_cula.d \
./src.cpp/carma_cusparse.d \
./src.cpp/carma_fft.d \
./src.cpp/carma_fft_conv.d \
./src.cpp/carma_host_obj.d \
./src.cpp/carma_ipcs.d \
./src.cpp/carma_ksparse.d \
./src.cpp/carma_magma.d \
./src.cpp/carma_multithread.d \
./src.cpp/carma_obj.d \
./src.cpp/carma_rng.d \
./src.cpp/carma_sparse_host_obj.d \
./src.cpp/carma_sparse_obj.d \
./src.cpp/carma_streams.d \
./src.cpp/carma_utils.d


# Each subdirectory must supply rules for building sources it contributes
src.cpp/%.o: ../src.cpp/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__GXX_EXPERIMENTAL_CXX0X__ -I/usr/local/cuda/include -O3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
