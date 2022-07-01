#pragma once

#include <string>

namespace mrpgo{;

class Posegraph;

class Loader{
public:
	bool Load    (std::string filename, Posegraph* pg);
	bool Save    (std::string filename, Posegraph* pg);
	void Sanitize(Posegraph* pg);

	Loader();
};

}
