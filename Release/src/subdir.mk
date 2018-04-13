################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BulkOnlySplines.cpp \
../src/BulkOnlySplinesQuadraticTail.cpp \
../src/ConfigItem.cpp \
../src/GaussianWavepacket.cpp \
../src/HeBulk.cpp \
../src/HeDrop.cpp \
../src/MathOperators.cpp \
../src/TDVMC.cpp \
../src/Utils.cpp 

OBJS += \
./src/BulkOnlySplines.o \
./src/BulkOnlySplinesQuadraticTail.o \
./src/ConfigItem.o \
./src/GaussianWavepacket.o \
./src/HeBulk.o \
./src/HeDrop.o \
./src/MathOperators.o \
./src/TDVMC.o \
./src/Utils.o 

CPP_DEPS += \
./src/BulkOnlySplines.d \
./src/BulkOnlySplinesQuadraticTail.d \
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


