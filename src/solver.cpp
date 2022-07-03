#include <solver.h>
#include <timer.h>

#include <limits>

namespace mrpgo{;

static Timer timer;
static Timer timer2;

static const double inf = std::numeric_limits<double>::max();
static const double eps = 1.0e-10;

///////////////////////////////////////////////////////////////////////////////////////////////////

Solver::Param::Param(){
	minStepSize     = 0.1;
	maxStepSize     = 1.0;
	cutoffStepSize  = 0.1;
	hastyStepSize   = false;
	regularization  = 0.001;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Solver::Status::Status(){
	obj       = inf;
	objDiff   = inf;
	stepSize  = 0.0;
	iterCount = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Solver::Solver(){
}

double Solver::CalcUpdatedObjective(double alpha){
	for(Variable* var : vars)
		var->Modify(alpha);

	return CalcObjective();
}

double Solver::CalcStepSize(){
	double amin = param.minStepSize;
	double amax = param.maxStepSize;
	double a[3], aprobe;
	double c[3], cprobe;
	
	if(amax - amin < param.cutoffStepSize)
		return 0.5 * (amin + amax);

	a[0] = amin;
	c[0] = CalcUpdatedObjective(amin);
	a[1] = a[2] = amax;
	c[1] = c[2] = CalcUpdatedObjective(amax);
	if(param.hastyStepSize && c[2] <= c[0])
		return a[2];

	if(c[2] > c[0]){
		while( c[1] > c[0] && a[1] - a[0] > param.cutoffStepSize * (amax - amin) ){
			a[1] = a[0] + 0.5 * (a[1] - a[0]);
			c[1] = CalcUpdatedObjective(a[1]);
		}
		if(param.hastyStepSize)
			return a[1];
	}

	while(a[2] - a[0] > param.cutoffStepSize){
		if(a[1] - a[0] > a[2] - a[1]){
			aprobe = 0.5 * (a[0] + a[1]);
			cprobe = CalcUpdatedObjective(aprobe);
			if(cprobe <= c[1]){
				a[2] = a[1];
				a[1] = aprobe;
				c[2] = c[1];
				c[1] = cprobe;
			}
			else{
				a[0] = aprobe;
				c[0] = cprobe;
			}
		}
		else{
			aprobe = 0.5 * (a[1] + a[2]);
			cprobe = CalcUpdatedObjective(aprobe);
			if(cprobe < c[1]){
				a[0] = a[1];
				a[1] = aprobe;
				c[0] = c[1];
				c[1] = cprobe;
			}
			else{
				a[2] = aprobe;
				c[2] = cprobe;
			}
		}
	}
	// returning average step size may cause drifting behavior in infeasible problems
	// returning lower value is safer
	//return 0.5 * (a[0] + a[2]);
	return a[0];
}

void Solver::Step(){
	for(auto& var : vars){
		var->ResetState();
	}

	timer.CountUS();
	CalcDirection();
	status.timeDir = timer.CountUS();

	timer.CountUS();
	status.stepSize = CalcStepSize();
	status.timeStep = timer.CountUS();

	for(auto& var : vars)
		var->Modify(status.stepSize);

   	double objPrev = status.obj;
	status.obj      = CalcObjective();
	status.objDiff  = status.obj - objPrev;
	status.iterCount++;

}

}
