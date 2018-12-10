#pragma once

#include <chrono>

using namespace std;

class Timer
{
private:
	bool isRunning;
	chrono::high_resolution_clock::time_point startTime;
	chrono::high_resolution_clock::time_point endTime;
	double lastIntervalTotalTime;

public:
	Timer()
	{
		isRunning = false;
		lastIntervalTotalTime = 0;
	}

	void start()
	{
		startTime = chrono::high_resolution_clock::now();
		isRunning = true;
	}

	void stop()
	{
		endTime = chrono::high_resolution_clock::now();
		isRunning = false;
	}

	double duration()
	{
		return chrono::duration_cast<chrono::milliseconds>((isRunning ? chrono::high_resolution_clock::now() : endTime) - startTime).count();
	}

	double durationSeconds()
	{
		return duration() / 1000.0;
	}

	double interval()
	{
		double tmp = duration();
		double interval = tmp - lastIntervalTotalTime;
		lastIntervalTotalTime = tmp;
		return interval;
	}
};
