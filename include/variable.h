#pragma once

#include <Eigen/Eigen>
using namespace Eigen;

namespace mrpgo{;

class Variable{
public:
	Vector3d  dx;            ///< delta
	
public:
	virtual void    ResetState(){
		dx = Vector3d(0.0, 0.0, 0.0);
	}
	virtual void	Modify(double rate) = 0;

	Variable(){}
};

template<class T>
class VariableImpl : public Variable{
public:
	T val;
	T val_tmp;

	virtual void ResetState(){
		Variable::ResetState();
		val_tmp = val;
	}

	VariableImpl():Variable(){}
};

/**
	scalar variable
 */
class SVar : public VariableImpl<double>{
public:
	virtual void Reset(){
		val = val_tmp = 0.0;
	}
	virtual void Modify(double alpha){
		val = val_tmp + alpha * dx[0];
	}

	SVar():VariableImpl(){
		Reset();
	}
};

/**
	2D vector variable
 */
class V2Var : public VariableImpl<Vector2d>{
public:
	virtual void Reset(){
		val     = Vector2d(0.0, 0.0);
		val_tmp = Vector2d(0.0, 0.0);
	}
	virtual void Modify(double alpha){
		val[0] = val_tmp[0] + alpha * dx[0];
		val[1] = val_tmp[1] + alpha * dx[1];
	}

	V2Var():VariableImpl(){
		Reset();
	}
};

/**
	3D vector variable
 */
class V3Var : public VariableImpl<Vector3d>{
public:
	virtual void Reset(){
		val     = Vector3d(0.0, 0.0, 0.0);
		val_tmp = Vector3d(0.0, 0.0, 0.0);
	}
	virtual void Modify(double alpha){
		val = val_tmp + alpha * dx;
	}

	V3Var():VariableImpl(){
		Reset();
	}
};

/**
	quaternion variable
 */
class QVar : public VariableImpl<Quaterniond>{
public:
	virtual void Reset(){
		val     = Quaterniond(0.0, 0.0, 0.0, 1.0);
		val_tmp = Quaterniond(0.0, 0.0, 0.0, 1.0);
	}
	virtual void Modify(double alpha){
		double dx_norm = dx.norm();
		val = AngleAxisd(alpha*dx_norm, dx/dx_norm) * val_tmp;
	}

	QVar():VariableImpl(){
		Reset();
	}
};

}
