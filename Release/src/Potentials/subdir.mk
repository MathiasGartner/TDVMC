################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Potentials/HFDB.cpp \
../src/Potentials/HFDB_He_He.cpp \
../src/Potentials/KTTY.cpp \
../src/Potentials/KTTY_He_Cs.cpp \
../src/Potentials/KTTY_He_Na.cpp \
../src/Potentials/LJ.cpp \
../src/Potentials/LJ_He_He.cpp \
../src/Potentials/LJ_XSplit.cpp \
../src/Potentials/SquareWell.cpp 

OBJS += \
./src/Potentials/HFDB.o \
./src/Potentials/HFDB_He_He.o \
./src/Potentials/KTTY.o \
./src/Potentials/KTTY_He_Cs.o \
./src/Potentials/KTTY_He_Na.o \
./src/Potentials/LJ.o \
./src/Potentials/LJ_He_He.o \
./src/Potentials/LJ_XSplit.o \
./src/Potentials/SquareWell.o 

CPP_DEPS += \
./src/Potentials/HFDB.d \
./src/Potentials/HFDB_He_He.d \
./src/Potentials/KTTY.d \
./src/Potentials/KTTY_He_Cs.d \
./src/Potentials/KTTY_He_Na.d \
./src/Potentials/LJ.d \
./src/Potentials/LJ_He_He.d \
./src/Potentials/LJ_XSplit.d \
./src/Potentials/SquareWell.d 


# Each subdirectory must supply rules for building sources it contributes
src/Potentials/%.o: ../src/Potentials/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++17 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


