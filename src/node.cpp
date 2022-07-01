#include <node.h>
#include <link.h>

namespace mrpgo{;

Node::Node(Posegraph* _owner){
	owner = _owner;
	index = -1;
}

void Node::CreateVariable(){
	int n  = owner->dim;
	int nv = owner->dim_trn;
	int nw = owner->dim_rot;
	
	A.Allocate(nv, nw);
	//b.Allocate(n);
	
	ddx   .Allocate(n);
	dx    .Allocate(n);
	dx_abs.Allocate(n);

	vec_clear(dx);
}

void Node::CalcJacobian(){
	Node* n0 = linkUpper->node[0];
	Node* n1 = linkUpper->node[1];

	if(owner->space == Posegraph::Space::SE2){
		r_upper.x() = n1->p.val.x() - n0->p.val.x();
		r_upper.y() = n1->p.val.y() - n0->p.val.y();
	}
	else{
		r_upper = n1->p.val - n0->p.val;
	}

	// flip sign if link is opposite
	if(n1 == this){
		//mat_vec_mul(linkUpper->ex_inv[1], linkUpper->e, b, -owner->biasFactor, 0.0);
	}
	else{
		//mat_vec_mul(linkUpper->ex_inv[0], linkUpper->e, b, -owner->biasFactor, 0.0);
		r_upper = -r_upper;
	}

	nodeUpper = linkUpper->Opposite(this);
    
	int rep = level + (depthBlock == 0 ? 0 : 1);
	while(rep--){
        Node* n1 = nodeUpper;
    
		/*
		vec_add(n1->b, b);
		if(owner->space == Posegraph::Space::SE2){
			double w = n1->b(2);
			b(0) += -r_upper.y()*w;
			b(1) +=  r_upper.x()*w;
		}
		else{
			vec3_t w;
			w.x = n1->b(3);
			w.y = n1->b(4);
			w.z = n1->b(5);
			vec3_t w_cross_r = w % r_upper;
			b(0) += w_cross_r.x();
			b(1) += w_cross_r.y();
			b(2) += w_cross_r.z();
		}
		*/

		r_upper += n1->r_upper;

		nodeUpper = nodeUpper->nodeUpper;
    }

	// calc A
    if(owner->space == Posegraph::Space::SE2){
		A(0,0) = -r_upper.y();
		A(1,0) =  r_upper.x();
	}
	else{
		Matrix3d rc;
		cross3_mat(-r_upper, rc);
		mat_copy(rc, A);
	}
}

}
