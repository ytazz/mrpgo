#pragma once

#include <corrector.h>

namespace mrpgo{;

class Link{
public:
	Node*        node[2];
	Posegraph*   owner;
		
	Vector3d     p_diff_abs;
	Vector3d     p_diff;
	Vector3d     p_ref;
	double       theta_diff;
	double       theta_ref;
	Quaterniond  q_diff;
	Quaterniond  q_ref;
	Matrix3d     R_ref;
		
	double     L;
	Matrix     ex[2];
	Matrix     ex_inv[2];
	Vector     Lx[2];
	Matrix     Lxx[2][2];
	Matrix     tmp;

	Matrix     I;
	Vector     e;
	Vector     d;

	Node* Opposite(Node* n);
	void  CreateVariable();
	void  CalcCost2     ();
	void  CalcCost3     ();
	void  CalcCostCommon();
	void  CalcCoef2     ();
	void  CalcCoef3     ();
	void  CalcCoefCommon();
		
	Link(Posegraph* _owner);
};

}
