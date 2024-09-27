#pragma once

#include <vector>

#include <blas.h>
#include <solver.h>
#include <posegraph.h>

#include <Eigen/Eigen>
using namespace Eigen;

namespace mrpgo{;

class Corrector : public Solver{
public:
	struct SolverType{
		enum{
			Custom,
			Cholmod,
		};
	};

	int       levelBase;
	int       minLevel;
	int       maxLevel;
	int       blockSize;
    double    maxRotation;
    double    regularization;
	int       numIter;
	double    shrinkRate;
	double    biasDecayRate;
	double    biasFactor;
	int       solverType;
	double    errorThreshold;
	int       numThreads;   ///< number of parallel threads.  0: use maximum available
	bool      verbose;
    
    double    J;
	double    maxError;

	int       Ttotal;

    struct Block{
		struct LxMap{
			Node*  node;
			Vector Lx;
			Vector Lx_upper;
		};
		struct LxxLxbMap{
			Node*  node;
			Matrix Lxx;
			Vector Lx;
		};
		struct LxxRowMap{
			Node*  node;
			Matrix Lxx;
			Matrix Lxx_upper;
		};
		struct LxxColMap{
			Node*  node;
			Matrix Lxx;
			Matrix Lxx_upper;
		};
		struct LxxLxMap{
			Node*  node;
			Vector Lx;
			Matrix Lxx;
		};
		struct LxBlockMap{
			Vector Lx;
			Vector Lx_block;
		};
		struct LxxBlockMap{
			Matrix Lxx;
			Matrix Lxx_block;
		};
		struct ddxBlockMap{
			Node*  node;
			Vector ddx_block;
		};
        Nodes             nodes;
        std::vector<int>  depths;
	    SparseMatrix      Lxx;   //< submatrix of Lxx corresponding to this block
	    SparseVector      Lx;    //< subvector of Lx
		SparseMatrix      Lxx_tmp;
		SparseVector      Lx_tmp;
		SparseVector      ddx;   //< subvector of ddx

		std::vector<LxMap>       Lx_map;
		std::vector<LxxLxbMap>   Lxx_Lx_b_map;
		std::vector<LxxRowMap>   Lxx_row_map;
		std::vector<LxxColMap>   Lxx_col_map;
		std::vector<LxxLxMap>    Lxx_Lx_map;
		std::vector<LxBlockMap> 	Lx_block_map;
		std::vector<LxxBlockMap> Lxx_block_map;
		std::vector<ddxBlockMap> ddx_block_map;

		LinearSolver*  linsolver;
	};
	struct Level{
        int                 m;
        std::vector<Block>  blocks;
    };

	Posegraph*          pg;
	bool                first;
	std::vector<Level>  levels;

	SparseVector    Lx ;
	SparseMatrix    Lxx;
	LinearSolver*   linsolver;

	std::string  statFilename;
	FILE*        fileStat;
	
public:
    void  SetSolverType (std::string str);
	void  Analyze       ();
	void  Correct       ();
	void  CalcCost      ();
	void  PrintSparsity (const char* filename, bool ordered);
	void  PrintDot      (const char* filename);

    virtual double  CalcObjective();
	virtual void    CalcDirection();

     Corrector(Posegraph* _pg);
	~Corrector();
};

}
