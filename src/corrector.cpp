#include <corrector.h>
#include <posegraph.h>
#include <node.h>
#include <link.h>
#include <timer.h>

#include <fstream>
#include <deque>

#include <omp.h>

namespace mrpgo{;

static Timer timer;
static Timer timer2;

using namespace std;

Corrector::Corrector(Posegraph* _pg){
	pg    = _pg;
	first = true;

	levelBase = 2;
	minLevel  = 0;
	maxLevel  = 0;
	blockSize = 100;
	
	maxRotation     = 1.0;
    regularization  = 0.001;
	numIter         = 1;
	shrinkRate      = 0.0;
	biasFactor      = 1.0;
	biasDecayRate   = 1.0;
	solverType      = SolverType::Custom;
	errorThreshold  = 0.0;
	numThreads      = 0;

	fileStat = 0;
}

Corrector::~Corrector(){
	if(fileStat)
		fclose(fileStat);
}

void Corrector::SetSolverType(string str){
	if(str == "custom" ) solverType = SolverType::Custom;
	if(str == "cholmod") solverType = SolverType::Cholmod;
}

void Corrector::Analyze(){
	deque<int> dtmp, dtmp_next;
	for(int d = 0; d < pg->depths.size(); d++)
		dtmp.push_back(d);

	levels.push_back(Level());
	int i = 0;
	while(!dtmp.empty()){
		Level& lv = levels[i];

		// check if sum is smaller than limit
		//int nsum = 0;
		//for(int d : dtmp)
		//	nsum += depths[d].nodes.size();

		// if smaller than limit or level reached maximum
		if(/*nsum < blockSize || */i == maxLevel){
			lv.blocks.push_back(Block());
			Block& bl = lv.blocks[0];
			for(int d : dtmp){
				for(Node* n : pg->depths[d].nodes){
					n->level = i;
					n->block = 0;
					n->indexBlock = bl.nodes .size();
					n->depthBlock = bl.depths.size();
					bl.nodes.push_back(n);
				}
				bl.depths.push_back(d);
			}
			break;
		}
		
		bool blockOpen = false;
		while(!dtmp.empty()){
			int d = dtmp.front();
			dtmp.pop_front();
			if(!blockOpen){
				dtmp_next.push_back(d);
				lv.blocks.push_back(Block());
				blockOpen = true;
				//nsum      = 0;
			}
			else{
				Block& bl = lv.blocks.back();
				for(Node* n : pg->depths[d].nodes){
					n->level = i;
					n->block = lv.blocks.size()-1;
					n->indexBlock = bl.nodes .size();
					n->depthBlock = bl.depths.size();
					bl.nodes.push_back(n);
				}
				bl.depths.push_back(d);
				//nsum += depths[d].nodes.size();
				if(/*nsum > blockSize || */bl.depths.size() == levelBase-1){
					blockOpen = false;
				}
			}
		}
		
		i++;
		dtmp.swap(dtmp_next);
		levels.push_back(Level());
	}

	maxLevel = levels.size()-1;
	
	// upper nodes in the hierarchy
	for(int i = 0; i < maxLevel; i++){
		Level& lv = levels[i];
		for(Block& bl : lv.blocks){
			for(int d : bl.depths){
				Posegraph::Depth& dp = pg->depths[d];
				for(Node* n : dp.nodes){
					n->nodeUpper = n->linkUpper->Opposite(n);
					
					int rep = n->level + (n->depthBlock == 0 ? 0 : 1);
					while(rep--){
						n->nodeUpper = n->nodeUpper->nodeUpper;
					}
    			}
			}
		}
	}
    
    for(int i = 0; i < levels.size(); i++){
        Level& lv = levels[i];
        //DSTR << "level " << i << endl;

        for(int j = 0; j < lv.blocks.size(); j++){
            Block& bl = lv.blocks[j];
            //DSTR << " block " << j << endl;
			//DSTR << "  depth ";
			//for(int d : bl.depths){
            //    DSTR << " " << d;
            //}
			//DSTR << endl;
        }
    }

	int N = pg->nodes.size();
   
	Lx .Resize(N, pg->dim);
	Lxx.Resize(N, N, pg->dim);

    for(Level& lv : levels){
        for(Block& bl : lv.blocks){
            int nx = bl.nodes.size();
            bl.Lxx.Resize(nx, nx, pg->dim);
		    bl.Lx .Resize(nx,     pg->dim);
		    bl.ddx.Resize(nx,     pg->dim);
        }
    }
    
	for(Node* n : pg->nodes)
		n->CreateVariable();

	for(Link* l : pg->links)
		l->CreateVariable();

	// max initial error
	CalcObjective();
	maxError = 0.0;
	for(Link* l : pg->links){
		maxError = std::max(maxError, l->L);
	}

	// initial guess
	for(int d = 1; d < pg->depths.size(); d++){
		for(Node* n1 : pg->depths[d].nodes){
			Link* l  = n1->linkUpper;
			Node* n0 = l->Opposite(n1);
			if(pg->space == Posegraph::Space::SE2){
				if(l->node[1] == n1){
					Vector2d tmp = n0->R2*Vector2d(l->p_ref.x(), l->p_ref.y());
					n1->p.val.x() = n0->p.val.x() + tmp.x();
					n1->p.val.y() = n0->p.val.y() + tmp.y();
					n1->theta.val = n0->theta.val + l->theta_ref;
					n1->R2 = Rotation2D<double>(n1->theta.val);
				}
				else{
					n1->theta.val = n0->theta.val - l->theta_ref;
					n1->R2 = Rotation2D<double>(n1->theta.val);
					Vector2d tmp = n1->R2*Vector2d(l->p_ref.x(), l->p_ref.y());
					n1->p.val.x() = n0->p.val.x() - tmp.x();
					n1->p.val.y() = n0->p.val.y() - tmp.y();
				}
			}
			else{
				if(l->node[1] == n1){
					n1->p.val = n0->p.val + n0->q.val*l->p_ref;
					n1->q.val = n0->q.val*l->q_ref;
				}
				else{
					n1->q.val = n0->q.val*l->q_ref.conjugate();
					n1->p.val = n0->p.val - n1->q.val*l->p_ref;
				}
				n1->R = n1->q.val.toRotationMatrix();
			}
		}
	}

	// allocate Lx and Lxx and store pointers
	for(Link* l : pg->links){
		int ix0  = l->node[0]->index;
		int ix1  = l->node[1]->index;
        
		l->Lx[0]     = Lx.SubVector(ix0);
		l->Lx[1]     = Lx.SubVector(ix1);
		l->Lxx[0][0] = Lxx.SubMatrix(ix0, ix0);
		l->Lxx[0][1] = Lxx.SubMatrix(ix0, ix1);
		l->Lxx[1][0] = Lxx.SubMatrix(ix1, ix0);
		l->Lxx[1][1] = Lxx.SubMatrix(ix1, ix1);
	}
	for(Node* n : pg->nodes){
		int    ix = n->index;
		Level& lv = levels[n->level];
		Block& bl = lv.blocks[n->block];
	
		n->Lxx = Lxx.SubMatrix(ix, ix);
	
		Block::ddxBlockMap m;
		m.node = n;
		m.ddx_block = bl.ddx.SubVector(n->indexBlock);
		bl.ddx_block_map.push_back(m);
	}
	
	for(int i = 0; i < maxLevel; i++){
		Level& lv = levels[i];
	
		for(int ix0 = 0; ix0 < Lxx.m; ix0++)for(auto it = Lxx.rows[ix0].begin(); it != Lxx.rows[ix0].end(); it++){
			int     ix1 = it->first ;
			Matrix& _m  = it->second;
			Node*   n0  = pg->nodes[ix0];
			Node*   n1  = pg->nodes[ix1];
	
			if(n1->level == i){
				Block&  bl = lv.blocks[n1->block];
				Block::LxxLxbMap m;
				m.node      = n1;
				m.Lxx       = _m;
				m.Lx        = Lx.SubVector(ix0);
				bl.Lxx_Lx_b_map.push_back(m);
			}
		}
		for(auto it = Lx.begin(); it != Lx.end(); it++){
			Node*   n  = pg->nodes[it->first];
	
			if(n->level == i){
				Block&  bl = lv.blocks[n->block];
				Block::LxMap m;
				m.node     = n;
				m.Lx       = it->second;
				m.Lx_upper = Lx.SubVector(n->nodeUpper->index);
				bl.Lx_map.push_back(m);
			}
		}

		for(int ix0 = 0; ix0 < Lxx.m; ix0++)for(auto it = Lxx.rows[ix0].begin(); it != Lxx.rows[ix0].end(); it++){
			int     ix1 = it->first ;
			Matrix& _m  = it->second;
			Node*   n0  = pg->nodes[ix0];
			Node*   n1  = pg->nodes[ix1];
	
			if(n0->level == i){
				Block&  bl = lv.blocks[n0->block];
				Block::LxxRowMap m;
				m.node      = n0;
				m.Lxx       = _m;
				m.Lxx_upper = Lxx.SubMatrix(n0->nodeUpper->index, ix1);
				bl.Lxx_row_map.push_back(m);
			}
		}
		for(int ix0 = 0; ix0 < Lxx.m; ix0++)for(auto it = Lxx.rows[ix0].begin(); it != Lxx.rows[ix0].end(); it++){
			int     ix1 = it->first ;
			Matrix& _m  = it->second;
			Node*   n0  = pg->nodes[ix0];
			Node*   n1  = pg->nodes[ix1];
	
			if(n1->level == i){
				Block&  bl = lv.blocks[n1->block];
				Block::LxxColMap m;
				m.node      = n1;
				m.Lxx       = _m;
				m.Lxx_upper = Lxx.SubMatrix(ix0, n1->nodeUpper->index);
				bl.Lxx_col_map.push_back(m);
			}
		}
	}
	for(int ix0 = 0; ix0 < Lxx.m; ix0++)for(auto it = Lxx.rows[ix0].begin(); it != Lxx.rows[ix0].end(); it++){
		int     ix1 = it->first ;
		Matrix& _m  = it->second;
		Node*   n0  = pg->nodes[ix0];
		Node*   n1  = pg->nodes[ix1];
		Level&  lv  = levels[n1->level];
		Block&  bl  = lv.blocks[n1->block];
	
		Block::LxxLxMap m;
		m.node = n1;
		m.Lxx  = _m;
		m.Lx   = Lx.SubVector(ix0);
		bl.Lxx_Lx_map.push_back(m);
	}
	
	for(auto it = Lx.begin(); it != Lx.end(); it++){
		Node*   n  = pg->nodes[it->first];
		Level&  lv = levels[n->level];
		Block&  bl = lv.blocks[n->block];
	
		Block::LxBlockMap m;
		m.Lx = it->second;
		m.Lx_block = bl.Lx.SubVector(n->indexBlock);
	
		bl.Lx_block_map.push_back(m);
	}
	for(int ix0 = 0; ix0 < Lxx.m; ix0++)for(auto it = Lxx.rows[ix0].begin(); it != Lxx.rows[ix0].end(); it++){
		int     ix1 = it->first ;
		Matrix& _m  = it->second;
		Node*   n0  = pg->nodes[ix0];
		Node*   n1  = pg->nodes[ix1];
	
		if( n0->level == n1->level && n0->block == n1->block){
			Level&  lv = levels[n0->level];
			Block&  bl = lv.blocks[n0->block];
	
			Block::LxxBlockMap m;
			m.Lxx = _m;
			m.Lxx_block = bl.Lxx.SubMatrix(n0->indexBlock, n1->indexBlock);
	
			bl.Lxx_block_map.push_back(m);
		}
	}

	// create linear solver
	for(Level& lv : levels){
		for(Block& bl : lv.blocks){
			if(solverType == SolverType::Custom)
				bl.linsolver = new LinearSolverCustom();
			if(solverType == SolverType::Cholmod)
				bl.linsolver = new LinearSolverCholmod();

			bl.linsolver->Init(bl.Lxx);
		}
	}

	vars.clear();
	if(pg->space == Posegraph::Space::SE2){
		for(Node* n : pg->nodes){
			vars.push_back(&n->p);
			vars.push_back(&n->theta);
		}
	}
	else{
		for(Node* n : pg->nodes){
			vars.push_back(&n->p);
			vars.push_back(&n->q);
		}
	}

	biasFactor = 1.0;
}

void Corrector::Correct(){
	if(first){
		// set number of parallel threads
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
			printf("using %d threads\n", numThreads);
		}
		else{
			printf("using maximum available threads: %d\n", omp_get_max_threads());
		}

		printf("analyzing\n");
		
		timer.CountUS();
		pg->Analyze();
		Analyze();
		int Tanalyze = timer.CountUS();

		printf(" elapsed time: %d [us]\n");
	
		if(!statFilename.empty()){
			fileStat = fopen(statFilename.c_str(), "w");
			if(fileStat){
				fprintf(fileStat, "Tpre1, Tpre2, Tpre3");
				for(int i = maxLevel; i >= 0; i--){
					fprintf(fileStat, ", Tsolve%d, Tupdate%d", i, i);
				}
				fprintf(fileStat, ", cost, step\n");
			}
		}
	
		first = false;
	}
	
	Solver::Step();
	
	if(fileStat)
		fprintf(fileStat, ", %f, %f\n", status.obj, status.stepSize);
	
	biasFactor *= biasDecayRate;
}

inline void cross2mattr_vec_mul_add(const Vector3d& r, const Vector& v, Vector& y){
	for(int i = 0; i < 3; i++)
		y(i) += v(i);
	y(2) += r.y()*v(0) - r.x()*v(1);
}
inline void cross2mattr_mat_mul_add(const Vector3d& r, const Matrix& m, Matrix& y){
	for(int i = 0; i < 3; i++)for(int j = 0; j < 3; j++)
		y(i,j) += m(i,j);
	for(int j = 0; j < 3; j++)
		y(2,j) += r.y()*m(0,j) - r.x()*m(1,j);
}
inline void mat_cross2mat_mul_add(const Vector3d& r, const Matrix& m, Matrix& y){
	for(int i = 0; i < 3; i++)for(int j = 0; j < 3; j++)
		y(i,j) += m(i,j);
	for(int i = 0; i < 3; i++)
		y(i,2) += m(i,0)*r.y() - m(i,1)*r.x();
}
inline void cross3mattr_vec_mul_add(const Vector3d& r, const Vector& v, Vector& y){
	for(int i = 0; i < 6; i++)
		y(i) += v(i);
	y(3) +=  r.z()*v(1) - r.y()*v(2);
	y(4) += -r.z()*v(0) + r.x()*v(2);
	y(5) +=  r.y()*v(0) - r.x()*v(1);
}
inline void cross3mattr_mat_mul_add(const Vector3d& r, const Matrix& m, Matrix& y){
	for(int i = 0; i < 6; i++)for(int j = 0; j < 6; j++)
		y(i,j) += m(i,j);
	for(int j = 0; j < 6; j++){
		y(3,j) +=  r.z()*m(1,j) - r.y()*m(2,j);
		y(4,j) += -r.z()*m(0,j) + r.x()*m(2,j);
		y(5,j) +=  r.y()*m(0,j) - r.x()*m(1,j);
	}
}
inline void mat_cross3mat_mul_add(const Vector3d& r, const Matrix& m, Matrix& y){
	for(int i = 0; i < 6; i++)for(int j = 0; j < 6; j++)
		y(i,j) += m(i,j);
	for(int i = 0; i < 6; i++){
		y(i,3) +=  r.z()*m(i,1) - r.y()*m(i,2);
		y(i,4) += -r.z()*m(i,0) + r.x()*m(i,2);
		y(i,5) +=  r.y()*m(i,0) - r.x()*m(i,1);
	}
}

void Corrector::CalcCost(){
    int tPrepare1, tPrepare2, tPrepare3;

    timer.CountUS();

	// clear gradient and hessian
	vec_clear(Lx );
    mat_clear(Lxx);
    
	// calc gradient and hessian of links
	if(pg->space == Posegraph::Space::SE2){
		for(auto& ld : pg->linksDisjoint){
			#pragma omp parallel for
			for(int i = 0; i < ld.size(); i++){
				Link* l = ld[i];
				l->CalcCost2();
				l->CalcCoef2();
				l->CalcCoefCommon();
			}
		}
	}
	else{
		for(auto& ld : pg->linksDisjoint){
			#pragma omp parallel for
			for(int i = 0; i < ld.size(); i++){
				Link* l = ld[i];
				l->CalcCost3();
				l->CalcCoef3();
				l->CalcCoefCommon();
			}
		}
	}

    tPrepare1 = timer.CountUS();
    //DSTR << "prepare1: " << tPrepare1 << endl;
    
    timer.CountUS();
    
    for(int i = 0; i < maxLevel; i++){
		Level& lv = levels[i];
		for(Block& bl : lv.blocks){
			for(int d : bl.depths){
				Posegraph::Depth& dp = pg->depths[d];
				
				#pragma omp parallel for
				for(int j = 0; j < dp.nodes.size(); j++){
					dp.nodes[j]->CalcJacobian();
				}
			}
		}
    }

    tPrepare2 = timer.CountUS();
    
    timer.CountUS();

	int nv = pg->dim_trn;
	int nw = pg->dim_rot;

    // setup lowest-level gradient and hessian
    {
        //#pragma omp parallel for
	    for(Node* n : pg->nodes){
            for(int i = 0; i < pg->dim; i++)
				n->Lxx(i,i) += regularization;
        }
    }

	// calc upper-level grad and hessian by aggrigation
	for(int i = 0; i < maxLevel; i++){
		Level& lv = levels[i];
	
#pragma omp parallel
		{
			if(pg->space == Posegraph::Space::SE2){
				#pragma omp for
				for(int ib = 0; ib < lv.blocks.size(); ib++){
					Block& bl = lv.blocks[ib];

					for(auto it = bl.Lx_map.begin(); it != bl.Lx_map.end(); it++){
						cross2mattr_vec_mul_add(-it->node->r_upper, it->Lx, it->Lx_upper);
					}
				}
		
				#pragma omp for
				for(int ib = 0; ib < lv.blocks.size(); ib++){
					Block& bl = lv.blocks[ib];

					for(auto it = bl.Lxx_row_map.begin(); it != bl.Lxx_row_map.end(); it++){
						cross2mattr_mat_mul_add(-it->node->r_upper, it->Lxx, it->Lxx_upper);
					}
				}
		
				#pragma omp for
				for(int ib = 0; ib < lv.blocks.size(); ib++){
					Block& bl = lv.blocks[ib];

					for(auto it = bl.Lxx_col_map.begin(); it != bl.Lxx_col_map.end(); it++){
						mat_cross2mat_mul_add(-it->node->r_upper, it->Lxx, it->Lxx_upper);
					}
				}
			}
			else{
				#pragma omp for
				for(int ib = 0; ib < lv.blocks.size(); ib++){
					Block& bl = lv.blocks[ib];

					for(auto it = bl.Lx_map.begin(); it != bl.Lx_map.end(); it++){
						cross3mattr_vec_mul_add(-it->node->r_upper, it->Lx, it->Lx_upper);
					}
				}
		
				#pragma omp for
				for(int ib = 0; ib < lv.blocks.size(); ib++){
					Block& bl = lv.blocks[ib];

					for(auto it = bl.Lxx_row_map.begin(); it != bl.Lxx_row_map.end(); it++){
						cross3mattr_mat_mul_add(-it->node->r_upper, it->Lxx, it->Lxx_upper);
					}
				}
		
				#pragma omp for
				for(int ib = 0; ib < lv.blocks.size(); ib++){
					Block& bl = lv.blocks[ib];

					for(auto it = bl.Lxx_col_map.begin(); it != bl.Lxx_col_map.end(); it++){
						mat_cross3mat_mul_add(-it->node->r_upper, it->Lxx, it->Lxx_upper);
					}
				}
			}
		}
    }

    tPrepare3 = timer.CountUS();
    
	if(fileStat)
		fprintf(fileStat, "%d, %d, %d", tPrepare1, tPrepare2, tPrepare3);
}

void Corrector::CalcDirection(){
	CalcCost();

	// cold start
	for(Node* n : pg->nodes){
		vec_clear(n->dx );	
		vec_clear(n->ddx);	
	}

	for(int n = 0; n < numIter; n++){
	    for(int i = maxLevel; i >= minLevel; i--){
            Level& lv  = levels[i];

			timer.CountUS();

			#pragma omp parallel for
            for(int j = 0; j < lv.blocks.size(); j++){
				Block& bl = lv.blocks[j];
                //timer2.CountUS();

		        // solve sparse linear equation of i-th level
				for(auto it = bl.Lx_block_map.begin(); it != bl.Lx_block_map.end(); it++){
					vec_copy(it->Lx, it->Lx_block);
				}
				for(auto it = bl.Lxx_block_map.begin(); it != bl.Lxx_block_map.end(); it++){
					mat_copy(it->Lxx, it->Lxx_block);
				}

				// only shallow copy
				bl.Lxx_tmp = bl.Lxx;
				bl.Lx_tmp  = bl.Lx;
				
				bl.linsolver->Solve(bl.Lxx_tmp, bl.Lx_tmp, bl.ddx);

		        // store result
				for(auto it = bl.ddx_block_map.begin(); it != bl.ddx_block_map.end(); it++){
					vec_copy(it->ddx_block, it->node->ddx);
					vec_add (it->ddx_block, it->node->dx );
				}
		    }
			int tSolve = timer.CountUS();

			timer.CountUS();
			for(int j = 0; j < lv.blocks.size(); j++){
				Block& bl = lv.blocks[j];
			
				for(auto it = bl.Lxx_Lx_map.begin(); it != bl.Lxx_Lx_map.end(); it++){
					mat_vec_mul(it->Lxx, it->node->ddx, it->Lx, 1.0, 1.0);
				}
			}
			int tUpdate = timer.CountUS();
			
			if(fileStat)
				fprintf(fileStat, ", %d, %d", tSolve, tUpdate);
	    }
    }
	
	int nv = pg->dim_trn;
	int nw = pg->dim_rot;

	// forward multi-res transform
	for(int i = maxLevel; i >= 0; i--){
        Level& lv = levels[i];
        for(int j = 0; j < lv.blocks.size(); j++){
            Block& bl = lv.blocks[j];
            for(Node* n : bl.nodes){
				int ix  = n->index;
				
				vec_copy(n->dx, n->dx_abs);

				if(i < maxLevel){
					//vec_add(n->b, n->dx_abs);
					vec_add(n->nodeUpper->dx_abs, n->dx_abs);
					mat_vec_mul(n->A, n->nodeUpper->dx_abs.SubVector(nv,nw), n->dx_abs.SubVector(0,nv), 1.0, 1.0);
				}
            }
        }
    }

	if(pg->space == Posegraph::Space::SE2){
		for(Node* n : pg->nodes){
			n->p    .dx[0] = n->dx_abs(0);
			n->p    .dx[1] = n->dx_abs(1);
			n->theta.dx[0] = n->dx_abs(2);
		}
	}
	else{
		for(Node* n : pg->nodes){
			n->p.dx[0] = n->dx_abs(0);
			n->p.dx[1] = n->dx_abs(1);
			n->p.dx[2] = n->dx_abs(2);
			n->q.dx[0] = n->dx_abs(3);
			n->q.dx[1] = n->dx_abs(4);
			n->q.dx[2] = n->dx_abs(5);
		}
	}
	
}

double Corrector::CalcObjective(){
    if(pg->space == Posegraph::Space::SE2){
		#pragma omp parallel for
		for(int i = 0; i < pg->nodes.size(); i++){
			Node* n = pg->nodes[i];
			n->R2 = Rotation2D<double>(n->theta.val);
		}
	}
	else{
		#pragma omp parallel for
		for(int i = 0; i < pg->nodes.size(); i++){
			Node* n = pg->nodes[i];
			n->R = n->q.val.toRotationMatrix();
		}
	}
    
	if(pg->space == Posegraph::Space::SE2){
		#pragma omp parallel for
		for(int i = 0; i < pg->links.size(); i++){
			Link* l = pg->links[i];
			l->CalcCost2();
		}
	}
	else{
		#pragma omp parallel for
		for(int i = 0; i < pg->links.size(); i++){
			Link* l = pg->links[i];
			l->CalcCost3();
		}
	}
	
	J = 0.0;
	for(Link* l : pg->links)
		J += l->L;

	//DSTR << "J " << J << endl;

	return J;
}


struct CustomPrintCallback : SparseMatrix::PrintCallback{
	Corrector*  owner;
	bool ordered;

	virtual bool Order(int lhs, int rhs){
		if(ordered){
			Node* n0 = owner->pg->nodes[lhs];
			Node* n1 = owner->pg->nodes[rhs];
			return (n0->level >  n1->level) ||
				   (n0->level == n1->level && n0->block <  n1->block) ||
				   (n0->level == n1->level && n0->block == n1->block && n0->index < n1->index);
		}
		return lhs < rhs;
	}
	virtual int  Label(int i, int j, bool nonzero){
		if(nonzero)
			return 1;

		if(owner->maxLevel != 0){
			Node* ni = owner->pg->nodes[i];
			Node* nj = owner->pg->nodes[j];
			if( ni->level == nj->level &&
				ni->block == nj->block )
				return ni->level + 2;  // 0 and 1 are for zero and nonzero
		}

		return 0;
	}
};

void Corrector::Print(bool ordered){
	int N = pg->nodes.size();

	CustomPrintCallback callback;
	callback.owner   = this;
	callback.ordered = ordered;
   
	ofstream file("hessian.dat");
	Lxx.PrintSparsity(file, &callback);
}

void Corrector::PrintDot(){
	ofstream file("posegraph.dot");

	file << "digraph Graph{" << endl;

	// nodes
	vector<string>  names(pg->nodes.size());

	stringstream ss;
	for(Node* n : pg->nodes){
		ss.str("");
		ss << "node" << n->index;
	
		names[n->index] = ss.str();
		file << ss.str() << " [label=\"\", shape=point, width=0.5, height=0.5]" << endl;
	}

	file << endl;

	// links
	for(Node* n : pg->nodes){
		if(!n->linkUpper)
			continue;
		Node* n0 = n->linkUpper->Opposite(n);

		file << names[n0->index] << " -> " << names[n->index] << " [color=black, penwidth=10.0, arrowhead=none]" << endl;
	}

	// non-tree links
	for(Link* l : pg->links){
		Node* n0 = l->node[0];
		Node* n1 = l->node[1];

		if(n0->linkUpper == l || n1->linkUpper == l)
			continue;
		
		file << names[n0->index] << " -> " << names[n1->index] << " [color=black, style=dashed, penwidth=2.0, arrowhead=none, constraint=false]" << endl;
	}

	file << "}" << endl;
}

}
