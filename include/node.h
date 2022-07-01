#pragma once

#include <posegraph.h>
#include <variable.h>
#include <blas.h>

namespace mrpgo{;

class Node{
public:
	Posegraph*      owner;
    Links           links;
    int             index;      //< index over all nodes
    int             indexBlock; //< index inside block
    int             depth;      //< index used for defining hierarchy
    int             block;
    int             level;      //< level that this node blongs to
	int             depthBlock; //< depth index in this level
	Link*           linkUpper;
    Node*           nodeUpper;
		
	//Vector3d     p;
	//Quaterniond  q;
	//Matrix3d     R;
	//double       theta;
	V3Var        p;
	QVar         q;
	SVar         theta;
	Matrix3d     R;
	Matrix2d     R2;
				     
	//Vector3d     dp;
	//Quaterniond  dq;
	//Vector2d     dp2;
	//double       dtheta;

	Vector3d  r_upper;
	Matrix    A;  ///< linear transform coefficients from upper nodes
	//Vector    b;
	
	Vector    dx;
	Vector    dx_abs;
	Vector    ddx;
	Matrix    Lxx;
        
	void  CreateVariable();
	void  CalcJacobian  ();
		
	Node(Posegraph* _owner);
};

}
