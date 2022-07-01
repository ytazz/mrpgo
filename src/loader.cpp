#include <loader.h>
#include <posegraph.h>
#include <node.h>
#include <link.h>

#include <csv.h>

using namespace std;

namespace mrpgo{;

Loader::Loader(){
	
}

bool Loader::Load(string filename, Posegraph* pg){
	CsvReader csv;
	if(!csv.Read(filename, " ", true)){
		printf("failed to load %s", filename.c_str());
		return false;
	}

	for(int i = 0; i < (int)csv.NumRow(); i++){
		if(csv.NumCol(i) < 1)
			continue;
		
		int j = 0;
		string tag = csv.Get<string>(i, j++);
		if(tag == "VERTEX_SE2"){
			// tag(1) id(1) pose(3)
			if(csv.NumCol(i) < 1+1+3)
				continue;
	
			Node* n = new Node(pg);
			n->index     = csv.Get<int>   (i, j++);
			n->p.val.x() = csv.Get<double>(i, j++);
			n->p.val.y() = csv.Get<double>(i, j++);
			n->theta.val = csv.Get<double>(i, j++);
			
			// recognize as 2D posegraph
			pg->space = Posegraph::Space::SE2;
			pg->nodes.push_back(n);
		}
		if(tag == "VERTEX_SE3:QUAT"){
			// tag(1) id(1) pose(7)
			if(csv.NumCol(i) < 1+1+7)
				continue;
			
			Node* n = new Node(pg);
			n->index     = csv.Get<int>   (i, j++);
			n->p.val.x() = csv.Get<double>(i, j++);
			n->p.val.y() = csv.Get<double>(i, j++);
			n->p.val.z() = csv.Get<double>(i, j++);
			n->q.val.x() = csv.Get<double>(i, j++);
			n->q.val.y() = csv.Get<double>(i, j++);
			n->q.val.z() = csv.Get<double>(i, j++);
			n->q.val.w() = csv.Get<double>(i, j++);
			
			// recognize as 3D posegraph
			pg->space = Posegraph::Space::SE3;
			pg->nodes.push_back(n);
		}
		if(tag == "EDGE_SE2"){
			// tag(1) id(2) pose(3) info(6)
			if(csv.NumCol(i) < 1+2+3+6)
				continue;

			int id0 = csv.Get<int>(i, j++);
			int id1 = csv.Get<int>(i, j++);

			Node* n0 = pg->nodes[id0];
			Node* n1 = pg->nodes[id1];
			
			Link* l = new Link(pg);
			l->node[0] = n0;
			l->node[1] = n1;
			pg->links.push_back(l);
			
			l->p_ref.x() = csv.Get<double>(i, j++);
			l->p_ref.y() = csv.Get<double>(i, j++);
			l->theta_ref = csv.Get<double>(i, j++);
			
			l->I.Allocate(3,3);
			for(int r = 0; r < 3; r++) for(int c = 0; c < 3; c++){
				l->I(r,c) = (r <= c ? csv.Get<double>(i, j++) : l->I(c,r));
			}
		}
		if(tag == "EDGE_SE3:QUAT"){
			// tag(1) id(2) pose(7) info(21)
			if(csv.NumCol(i) < 1+2+7+21)
				continue;
			
			int id0 = csv.Get<int>(i, j++);
			int id1 = csv.Get<int>(i, j++);

			Node* n0 = pg->nodes[id0];
			Node* n1 = pg->nodes[id1];
			
			Link* l = new Link(pg);
			l->node[0] = n0;
			l->node[1] = n1;
			pg->links.push_back(l);
			
			l->p_ref.x() = csv.Get<double>(i, j++);
			l->p_ref.y() = csv.Get<double>(i, j++);
			l->p_ref.z() = csv.Get<double>(i, j++);
			l->q_ref.x() = csv.Get<double>(i, j++);
			l->q_ref.y() = csv.Get<double>(i, j++);
			l->q_ref.z() = csv.Get<double>(i, j++);
			l->q_ref.w() = csv.Get<double>(i, j++);

			l->I.Allocate(6,6);
			for(int r = 0; r < 6; r++) for(int c = 0; c < 6; c++)
				l->I(r,c) = (r <= c ? csv.Get<double>(i, j++) : l->I(c,r));

		}
	}

	return true;
}

bool Loader::Save(string filename, Posegraph* pg){
	FILE* file = fopen(filename.c_str(), "w");

	for(int i = 0; i < pg->nodes.size(); i++){
		Node* n = pg->nodes[i];

		fprintf(file, "%s", pg->space == Posegraph::Space::SE2 ? "VERTEX_SE2" : "VERTEX_SE3:QUAT");
		fprintf(file, " %d", i);

		if(pg->space == Posegraph::Space::SE2){
			fprintf(file, " %6.6f %6.6f %6.6f",
				n->p.val.x(),
				n->p.val.y(),
				n->theta.val
			);
		}
		else{
			fprintf(file, " %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
				n->p.val.x(),
				n->p.val.y(),
				n->p.val.z(),
				n->q.val.x(),
				n->q.val.y(),
				n->q.val.z(),
				n->q.val.w()
			);
		}
		fprintf(file, "\n");
	}

	for(Link* l : pg->links){
		fprintf(file, "%s", pg->space == Posegraph::Space::SE2 ? "EDGE_SE2" : "EDGE_SE3:QUAT");
		fprintf(file, " %d %d", l->node[0]->index, l->node[1]->index);
		
		if(pg->space == Posegraph::Space::SE2){
			fprintf(file, " %6.6f %6.6f %6.6f",
				l->p_ref.x(),
				l->p_ref.y(),
				l->theta_ref
			);
		}
		else{
			fprintf(file, " %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
				l->p_ref.x(),
				l->p_ref.y(),
				l->p_ref.z(),
				l->q_ref.x(),
				l->q_ref.y(),
				l->q_ref.z(),
				l->q_ref.w()
			);
		}

		if(pg->space == Posegraph::Space::SE2){
			fprintf(file, " %6.6f %6.6f %6.6f"
				                " %6.6f %6.6f"
				                      " %6.6f",
				l->I(0,0), l->I(0,1), l->I(0,2),
				           l->I(1,1), l->I(1,2),
				                      l->I(2,2)
			);
		}
		else{
			fprintf(file, " %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f"
				                " %6.6f %6.6f %6.6f %6.6f %6.6f"
				                      " %6.6f %6.6f %6.6f %6.6f"
				                            " %6.6f %6.6f %6.6f"
				                                  " %6.6f %6.6f"
				                                        " %6.6f",
				l->I(0,0), l->I(0,1), l->I(0,2), l->I(0,3), l->I(0,4), l->I(0,5),
				           l->I(1,1), l->I(1,2), l->I(1,3), l->I(1,4), l->I(1,5),
				                      l->I(2,2), l->I(2,3), l->I(2,4), l->I(2,5),
				                                 l->I(3,3), l->I(3,4), l->I(3,5),
				                                            l->I(4,4), l->I(4,5),
				                                                       l->I(5,5)
			);
		}
		fprintf(file, "\n");
	}

	fclose(file);

	return true;
}

void Loader::Sanitize(Posegraph* pg){
	// check for non-positive definite hessians
	Vector wr(6);
	Vector wi(6);
	Matrix vl(6,6);
	Matrix vr(6,6);
	Matrix Wr(6,6);
	Matrix tmp(6,6);
    double wmin, wmax;

	Matrix Ibase(6,6);

	for(Link* l : pg->links){
		/*
		// check positive-definiteness of information matrix
		mat_eig(Ibase, wr, wi, vl,vr);
        wmin =  inf;
        wmax = -inf;
		for(int i = 0; i < 6; i++){
            wmin = std::min(wmin, wr.vh[i]);
            wmax = std::max(wmax, wr.vh[i]);
			if(wr.vh[i] < 0.0){
				DSTR << "negative eigenvalue!" << endl;
				wr.vh[i] = 0.0;
			}
		}
        //DSTR << "wmin: " << wmin << " wmax: " << wmax << endl;
		mat_diag(wr, Wr);
		mat_mat_mul  (vr, Wr, tmp, 1.0f, 0.0f);
		mat_mattr_mul(tmp, vr, Ibase, 1.0f, 0.0f);
		//l->Ibase = vr*diag(wr)*vr.trans();
		*/

		// assign constant hessian (temporary code!)
		mat_clear(l->I);
		for(int i = 0; i < 6; i++)
			l->I(i,i) = 100.0;

	}

}

}
