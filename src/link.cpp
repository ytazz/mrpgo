#include <link.h>
#include <node.h>

namespace mrpgo{;

const double pi = 3.1415926535;

template<class T>
inline T WrapRadian(T v){
	while(v > (T)pi)
		v -= (T)(2.0 * pi);
	while(v < -pi)
		v += (T)(2.0 * pi);
	return v;
}

Link::Link(Posegraph* _owner){
	owner      = _owner;
	node[0]    = 0;
	node[1]    = 0;
	theta_ref  = 0.0;
	L          = 0.0;
}

Node* Link::Opposite(Node* n){
    return node[0] == n ? node[1] : node[0];
}
		
void Link::CreateVariable(){
	int n = owner->dim;

    e        .Allocate(n);
	d        .Allocate(n);
	ex    [0].Allocate(n, n);
	ex    [1].Allocate(n, n);
	ex_inv[0].Allocate(n, n);
	ex_inv[1].Allocate(n, n);
    tmp      .Allocate(n, n);
}

void Link::CalcCost2(){
	Node* n0 = node[0];
	Node* n1 = node[1];

	p_diff_abs.x() = n1->p.val.x() - n0->p.val.x();
	p_diff_abs.y() = n1->p.val.y() - n0->p.val.y();
	Vector2d tmp = n0->R2.transpose()*Vector2d(p_diff_abs.x(), p_diff_abs.y());
	p_diff.x() = tmp.x();
	p_diff.y() = tmp.y();
	theta_diff = n1->theta.val - n0->theta.val;

	e(0) = (p_diff.x() - p_ref.x());
	e(1) = (p_diff.y() - p_ref.y());
	e(2) = WrapRadian(theta_diff - theta_ref);

	CalcCostCommon();
}

void Link::CalcCost3(){
	Node* n0 = node[0];
	Node* n1 = node[1];

	p_diff_abs = n1->p.val - n0->p.val;
	p_diff = n0->q.val.conjugate()*p_diff_abs;
	q_diff = n0->q.val.conjugate()* n1->q.val;
			
	Vector3d ev = p_diff - p_ref;
	//vec3_t ew = q_ref*ToAxisAngle(q_ref.Conjugated()*q_diff);
	AngleAxisd tmp(q_ref.conjugate()*q_diff);
	Vector3d ew = 0.5*tmp.angle()*tmp.axis();
	//quat_t q_error = q_ref.Conjugated()*q_diff;
	//vec3_t ew = q_ref*(q_error.x, q_error.y, q_error.z);
	e(0) = ev[0];
	e(1) = ev[1];
	e(2) = ev[2];
	e(3) = ew[0];
	e(4) = ew[1];
	e(5) = ew[2];

	CalcCostCommon();
}

void Link::CalcCostCommon(){		
	mat_vec_mul(I, e, d, 1.0, 0.0);
	//L = (1.0/2.0)*vec_dot(d, e);
	L = vec_dot(d, e);
}

void Link::CalcCoef2(){
	Node* n0 = node[0];
	
	ex[0](0,0) = (-n0->R2(0,0)); ex[0](0,1) = (-n0->R2(1,0)); ex[0](0,2) = ( p_diff.y());
	ex[0](1,0) = (-n0->R2(0,1)); ex[0](1,1) = (-n0->R2(1,1)); ex[0](1,2) = (-p_diff.x());
	ex[0](2,0) = ( 0.0        ); ex[0](2,1) = ( 0.0        ); ex[0](2,2) = (-1.0       );
	
	ex[1](0,0) = ( n0->R2(0,0)); ex[1](0,1) = ( n0->R2(1,0)); ex[1](0,2) = ( 0.0     );
	ex[1](1,0) = ( n0->R2(0,1)); ex[1](1,1) = ( n0->R2(1,1)); ex[1](1,2) = ( 0.0     );
	ex[1](2,0) = ( 0.0        ); ex[1](2,1) = ( 0.0        ); ex[1](2,2) = ( 1.0     );
	
	ex_inv[0](0,0) = (-n0->R2(0,0)); ex_inv[0](0,1) = (-n0->R2(0,1)); ex_inv[0](0,2) = (-p_diff_abs.y());
	ex_inv[0](1,0) = (-n0->R2(1,0)); ex_inv[0](1,1) = (-n0->R2(1,1)); ex_inv[0](1,2) = ( p_diff_abs.x());
	ex_inv[0](2,0) = ( 0.0        ); ex_inv[0](2,1) = ( 0.0        ); ex_inv[0](2,2) = (-1.0           );
	
	ex_inv[1](0,0) = ( n0->R2(0,0)); ex_inv[1](0,1) = ( n0->R2(0,1)); ex_inv[1](0,2) = ( 0.0     );
	ex_inv[1](1,0) = ( n0->R2(1,0)); ex_inv[1](1,1) = ( n0->R2(1,1)); ex_inv[1](1,2) = ( 0.0     );
	ex_inv[1](2,0) = ( 0.0        ); ex_inv[1](2,1) = ( 0.0        ); ex_inv[1](2,2) = ( 1.0     );
}

void Link::CalcCoef3(){
	Node* n0 = node[0];
	Node* n1 = node[1];

	Quaterniond qtmp = n0->q.val * q_ref;
	Matrix3d Rtmp = qtmp.toRotationMatrix();
	
	Matrix3d rc;
	cross3_mat(p_diff_abs, rc);
	
	mat_copy (Matrix3d(-n0->R.transpose())   , ex[0].SubMatrix(0,0,3,3));
	mat_copy (Matrix3d( n0->R.transpose()*rc), ex[0].SubMatrix(0,3,3,3));
	mat_clear(                                 ex[0].SubMatrix(3,0,3,3));
	mat_copy (Matrix3d(-0.5*Rtmp.transpose()), ex[0].SubMatrix(3,3,3,3));

	mat_copy (Matrix3d( n0->R.transpose()   ), ex[1].SubMatrix(0,0,3,3));
	mat_clear(                                 ex[1].SubMatrix(0,3,3,3));
	mat_clear(                                 ex[1].SubMatrix(3,0,3,3));
	mat_copy (Matrix3d( 0.5*Rtmp.transpose()), ex[1].SubMatrix(3,3,3,3));

	mat_copy (Matrix3d(-n0->R        ), ex_inv[0].SubMatrix(0,0,3,3));
	mat_copy (Matrix3d(-rc*(2.0*Rtmp)), ex_inv[0].SubMatrix(0,3,3,3));
	mat_clear(                          ex_inv[0].SubMatrix(3,0,3,3));
	mat_copy (Matrix3d(-2.0*Rtmp     ), ex_inv[0].SubMatrix(3,3,3,3));
	
	mat_copy (Matrix3d( n0->R   ), ex_inv[1].SubMatrix(0,0,3,3));
	mat_clear(                     ex_inv[1].SubMatrix(0,3,3,3));
	mat_clear(                     ex_inv[1].SubMatrix(3,0,3,3));
	mat_copy (Matrix3d( 2.0*Rtmp), ex_inv[1].SubMatrix(3,3,3,3));
}

void Link::CalcCoefCommon(){

	// add to global gradient and hessian
	mattr_vec_mul(ex[0], d, Lx[0], 1.0, 1.0);
	mattr_vec_mul(ex[1], d, Lx[1], 1.0, 1.0);
	mattr_mat_mul(ex[0], I, tmp, 1.0, 0.0);
	mat_mat_mul(tmp, ex[0], Lxx[0][0], 1.0, 1.0);
	mat_mat_mul(tmp, ex[1], Lxx[0][1], 1.0, 1.0);
	mattr_mat_mul(ex[1], I, tmp, 1.0, 0.0);
	mat_mat_mul(tmp, ex[0], Lxx[1][0], 1.0, 1.0);
	mat_mat_mul(tmp, ex[1], Lxx[1][1], 1.0, 1.0);
}

}
