#pragma once
#include <chrono>
#include <iostream>

class HighResStopwatch
{
public:
	HighResStopwatch();
	~HighResStopwatch();
	//this is not quite as exact as including timer hardcode, but is close enough for keeping library

	bool timerRunning;
	std::chrono::time_point<std::chrono::system_clock, std::chrono::system_clock::duration> start;
	std::chrono::time_point<std::chrono::system_clock, std::chrono::system_clock::duration> stop;

	void StartTimer();
	void StopTimer();
	void Report();
};

