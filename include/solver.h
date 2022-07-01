#pragma once

#include <variable.h>

namespace mrpgo{;

class Solver{
public:
	struct Param{
		bool     verbose;
		double   minStepSize;
		double   maxStepSize;
		double   cutoffStepSize;
		bool     hastyStepSize;
		double   regularization;
		
		Param();
	};

	struct Status{
		double  obj;        ///< value of objective function
		double  objDiff;    ///< change of value of objective function
		double  stepSize;   ///< step size
		int     iterCount;  ///< cumulative iteration count
		int     timeDir;
		int     timeStep;

		Status();
	};

	Param            param;
	Status           status;
	
	typedef std::vector<Variable*>  Variables;
	Variables        vars;			///< array of all variables
	
public:
	double  CalcUpdatedObjective(double alpha);
	double  CalcStepSize();
	void    Step();
	
	/// evaluate objective function
	virtual double  CalcObjective() = 0;

	/// calculate update direction
	virtual void    CalcDirection() = 0;
	
public:

	Solver();

};

}
