#pragma once

#include <chrono>

namespace mrpgo{;

class Timer{
public:
	std::chrono::high_resolution_clock::time_point  last;

	/// count time in [us]
	int CountUS();

	Timer();
	~Timer();
};

}
