################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Observables/Grid.cpp \
../src/Observables/Observable.cpp \
../src/Observables/ObservableCollection.cpp \
../src/Observables/ObservableV.cpp \
../src/Observables/ObservableVsOnGrid.cpp \
../src/Observables/ObservableVsOnGridWithScaling.cpp \
../src/Observables/ObservableVsOnMultiGrid.cpp 

OBJS += \
./src/Observables/Grid.o \
./src/Observables/Observable.o \
./src/Observables/ObservableCollection.o \
./src/Observables/ObservableV.o \
./src/Observables/ObservableVsOnGrid.o \
./src/Observables/ObservableVsOnGridWithScaling.o \
./src/Observables/ObservableVsOnMultiGrid.o 

CPP_DEPS += \
./src/Observables/Grid.d \
./src/Observables/Observable.d \
./src/Observables/ObservableCollection.d \
./src/Observables/ObservableV.d \
./src/Observables/ObservableVsOnGrid.d \
./src/Observables/ObservableVsOnGridWithScaling.d \
./src/Observables/ObservableVsOnMultiGrid.d 


# Each subdirectory must supply rules for building sources it contributes
src/Observables/%.o: ../src/Observables/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


