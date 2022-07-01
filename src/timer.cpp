#include <timer.h>

namespace mrpgo{;

Timer::Timer(){
}

Timer::~Timer(){
}

int Timer::CountUS(){
	auto t = std::chrono::high_resolution_clock::now();
	auto duration = t - last;
	int us = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();

	last = t;

	return us;
}

}
