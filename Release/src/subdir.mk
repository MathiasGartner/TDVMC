################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BulkOnlySplines.cpp \
../src/BulkOnlySplinesOriginal.cpp \
../src/BulkQT.cpp \
../src/BulkSplines.cpp \
../src/BulkSplinesInterface.cpp \
../src/BulkSplinesNamespace.cpp \
../src/ConfigItem.cpp \
../src/GaussianWavepacket.cpp \
../src/HeBulk.cpp \
../src/HeDrop.cpp \
../src/MathOperators.cpp \
../src/TDVMC.cpp \
../src/Utils.cpp 

OBJS += \
./src/BulkOnlySplines.o \
./src/BulkOnlySplinesOriginal.o \
./src/BulkQT.o \
./src/BulkSplines.o \
./src/BulkSplinesInterface.o \
./src/BulkSplinesNamespace.o \
./src/ConfigItem.o \
./src/GaussianWavepacket.o \
./src/HeBulk.o \
./src/HeDrop.o \
./src/MathOperators.o \
./src/TDVMC.o \
./src/Utils.o 

CPP_DEPS += \
./src/BulkOnlySplines.d \
./src/BulkOnlySplinesOriginal.d \
./src/BulkQT.d \
./src/BulkSplines.d \
./src/BulkSplinesInterface.d \
./src/BulkSplinesNamespace.d \
./src/ConfigItem.d \
./src/GaussianWavepacket.d \
./src/HeBulk.d \
./src/HeDrop.d \
./src/MathOperators.d \
./src/TDVMC.d \
./src/Utils.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -I/usr/include/jsoncpp -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


