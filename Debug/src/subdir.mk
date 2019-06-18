################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CSDataBulkSplines.cpp \
../src/ConfigItem.cpp \
../src/CorrelationFunctionData.cpp \
../src/MathOperators.cpp \
../src/OneParticleData.cpp \
../src/ParticlePairProperties.cpp \
../src/ParticleProperties.cpp \
../src/SimulationStepData.cpp \
../src/SplineFactory.cpp \
../src/TDVMC.cpp \
../src/Utils.cpp 

OBJS += \
./src/CSDataBulkSplines.o \
./src/ConfigItem.o \
./src/CorrelationFunctionData.o \
./src/MathOperators.o \
./src/OneParticleData.o \
./src/ParticlePairProperties.o \
./src/ParticleProperties.o \
./src/SimulationStepData.o \
./src/SplineFactory.o \
./src/TDVMC.o \
./src/Utils.o 

CPP_DEPS += \
./src/CSDataBulkSplines.d \
./src/ConfigItem.d \
./src/CorrelationFunctionData.d \
./src/MathOperators.d \
./src/OneParticleData.d \
./src/ParticlePairProperties.d \
./src/ParticleProperties.d \
./src/SimulationStepData.d \
./src/SplineFactory.d \
./src/TDVMC.d \
./src/Utils.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


