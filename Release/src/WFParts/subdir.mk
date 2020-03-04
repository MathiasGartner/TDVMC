################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/WFParts/PairCorrelation.cpp \
../src/WFParts/SingleParticleFunction.cpp \
../src/WFParts/SplinedFunction.cpp 

OBJS += \
./src/WFParts/PairCorrelation.o \
./src/WFParts/SingleParticleFunction.o \
./src/WFParts/SplinedFunction.o 

CPP_DEPS += \
./src/WFParts/PairCorrelation.d \
./src/WFParts/SingleParticleFunction.d \
./src/WFParts/SplinedFunction.d 


# Each subdirectory must supply rules for building sources it contributes
src/WFParts/%.o: ../src/WFParts/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++17 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


