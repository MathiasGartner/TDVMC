################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Potentials/HFDB_He_He.cpp \
../src/Potentials/KTTY.cpp \
../src/Potentials/KTTY_He_Cs.cpp \
../src/Potentials/KTTY_He_Na.cpp \
../src/Potentials/LJ_He_He.cpp 

OBJS += \
./src/Potentials/HFDB_He_He.o \
./src/Potentials/KTTY.o \
./src/Potentials/KTTY_He_Cs.o \
./src/Potentials/KTTY_He_Na.o \
./src/Potentials/LJ_He_He.o 

CPP_DEPS += \
./src/Potentials/HFDB_He_He.d \
./src/Potentials/KTTY.d \
./src/Potentials/KTTY_He_Cs.d \
./src/Potentials/KTTY_He_Na.d \
./src/Potentials/LJ_He_He.d 


# Each subdirectory must supply rules for building sources it contributes
src/Potentials/%.o: ../src/Potentials/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


