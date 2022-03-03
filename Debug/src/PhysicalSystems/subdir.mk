################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/PhysicalSystems/BosonCluster.cpp \
../src/PhysicalSystems/BosonClusterWithLog.cpp \
../src/PhysicalSystems/BosonClusterWithLogParam.cpp \
../src/PhysicalSystems/BosonMixtureCluster.cpp \
../src/PhysicalSystems/BosonMixtureCluster_4thorder.cpp \
../src/PhysicalSystems/Bosons1D.cpp \
../src/PhysicalSystems/Bosons1D0th.cpp \
../src/PhysicalSystems/Bosons1D4th.cpp \
../src/PhysicalSystems/Bosons1DMixture.cpp \
../src/PhysicalSystems/Bosons1DSp.cpp \
../src/PhysicalSystems/BosonsBulk.cpp \
../src/PhysicalSystems/BosonsBulkDamped.cpp \
../src/PhysicalSystems/BosonsDiscrete.cpp \
../src/PhysicalSystems/BulkOnlySplines.cpp \
../src/PhysicalSystems/BulkOnlySplinesOriginal.cpp \
../src/PhysicalSystems/BulkSplines.cpp \
../src/PhysicalSystems/BulkSplinesPhi.cpp \
../src/PhysicalSystems/BulkSplinesScaled.cpp \
../src/PhysicalSystems/GaussianWavepacket.cpp \
../src/PhysicalSystems/HeBulk.cpp \
../src/PhysicalSystems/HeDrop.cpp \
../src/PhysicalSystems/InhBosons1D.cpp \
../src/PhysicalSystems/InhContactBosons.cpp \
../src/PhysicalSystems/LinearChain.cpp \
../src/PhysicalSystems/NUBosonsBulk.cpp \
../src/PhysicalSystems/NUBosonsBulkPB.cpp \
../src/PhysicalSystems/NUBosonsBulkPBBox.cpp \
../src/PhysicalSystems/NUBosonsBulkPBBoxAndRadial.cpp \
../src/PhysicalSystems/NUBosonsBulkPBWhitehead.cpp 

OBJS += \
./src/PhysicalSystems/BosonCluster.o \
./src/PhysicalSystems/BosonClusterWithLog.o \
./src/PhysicalSystems/BosonClusterWithLogParam.o \
./src/PhysicalSystems/BosonMixtureCluster.o \
./src/PhysicalSystems/BosonMixtureCluster_4thorder.o \
./src/PhysicalSystems/Bosons1D.o \
./src/PhysicalSystems/Bosons1D0th.o \
./src/PhysicalSystems/Bosons1D4th.o \
./src/PhysicalSystems/Bosons1DMixture.o \
./src/PhysicalSystems/Bosons1DSp.o \
./src/PhysicalSystems/BosonsBulk.o \
./src/PhysicalSystems/BosonsBulkDamped.o \
./src/PhysicalSystems/BosonsDiscrete.o \
./src/PhysicalSystems/BulkOnlySplines.o \
./src/PhysicalSystems/BulkOnlySplinesOriginal.o \
./src/PhysicalSystems/BulkSplines.o \
./src/PhysicalSystems/BulkSplinesPhi.o \
./src/PhysicalSystems/BulkSplinesScaled.o \
./src/PhysicalSystems/GaussianWavepacket.o \
./src/PhysicalSystems/HeBulk.o \
./src/PhysicalSystems/HeDrop.o \
./src/PhysicalSystems/InhBosons1D.o \
./src/PhysicalSystems/InhContactBosons.o \
./src/PhysicalSystems/LinearChain.o \
./src/PhysicalSystems/NUBosonsBulk.o \
./src/PhysicalSystems/NUBosonsBulkPB.o \
./src/PhysicalSystems/NUBosonsBulkPBBox.o \
./src/PhysicalSystems/NUBosonsBulkPBBoxAndRadial.o \
./src/PhysicalSystems/NUBosonsBulkPBWhitehead.o 

CPP_DEPS += \
./src/PhysicalSystems/BosonCluster.d \
./src/PhysicalSystems/BosonClusterWithLog.d \
./src/PhysicalSystems/BosonClusterWithLogParam.d \
./src/PhysicalSystems/BosonMixtureCluster.d \
./src/PhysicalSystems/BosonMixtureCluster_4thorder.d \
./src/PhysicalSystems/Bosons1D.d \
./src/PhysicalSystems/Bosons1D0th.d \
./src/PhysicalSystems/Bosons1D4th.d \
./src/PhysicalSystems/Bosons1DMixture.d \
./src/PhysicalSystems/Bosons1DSp.d \
./src/PhysicalSystems/BosonsBulk.d \
./src/PhysicalSystems/BosonsBulkDamped.d \
./src/PhysicalSystems/BosonsDiscrete.d \
./src/PhysicalSystems/BulkOnlySplines.d \
./src/PhysicalSystems/BulkOnlySplinesOriginal.d \
./src/PhysicalSystems/BulkSplines.d \
./src/PhysicalSystems/BulkSplinesPhi.d \
./src/PhysicalSystems/BulkSplinesScaled.d \
./src/PhysicalSystems/GaussianWavepacket.d \
./src/PhysicalSystems/HeBulk.d \
./src/PhysicalSystems/HeDrop.d \
./src/PhysicalSystems/InhBosons1D.d \
./src/PhysicalSystems/InhContactBosons.d \
./src/PhysicalSystems/LinearChain.d \
./src/PhysicalSystems/NUBosonsBulk.d \
./src/PhysicalSystems/NUBosonsBulkPB.d \
./src/PhysicalSystems/NUBosonsBulkPBBox.d \
./src/PhysicalSystems/NUBosonsBulkPBBoxAndRadial.d \
./src/PhysicalSystems/NUBosonsBulkPBWhitehead.d 


# Each subdirectory must supply rules for building sources it contributes
src/PhysicalSystems/%.o: ../src/PhysicalSystems/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++17 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


