#include "HighResStopwatch.h"


HighResStopwatch::HighResStopwatch()
{
	timerRunning = false;
	start = std::chrono::high_resolution_clock::now();
	stop = std::chrono::high_resolution_clock::now();
}


HighResStopwatch::~HighResStopwatch()
{
}

void HighResStopwatch::StartTimer()
{
	if (timerRunning)
	{
		std::cout << "Timer had already started" << std::endl << "This will delay original timer" << std::endl;
	}
	else
	{
		timerRunning = true;
		start = std::chrono::high_resolution_clock::now();
	}

}

void HighResStopwatch::StopTimer()
{
	stop = std::chrono::high_resolution_clock::now();

	if (!timerRunning)
	{
		std::cout << "Timer was not running" << std::endl
			<< "This will have effected any previous timer result" << std::endl;
	}
	else
	{
		timerRunning = false;
	}
}

void HighResStopwatch::Report()
{
	std::cout << "Timer recorded duration of "
		<< std::chrono::duration<double>(stop - start).count()
		<< " seconds" << std::endl;
}