#include <blas.h>

#include <algorithm>

namespace mrpgo{

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////

Vector::Vector(){
	n      = 0;
	vh     = 0;
}

Vector::Vector(int _n){
	n      = 0;
	vh     = 0;
	Allocate(_n);
}

Vector::~Vector(){
	Delete();
}

Vector Vector::SubVector(int ofst, int sz){
	Vector y;
	y.vh = vh + ofst;
	y.n  = sz;
	return y;
}

void Vector::Delete(){
	//if(vh) delete[] vh;
}

void Vector::Allocate(int _n){
	Delete();
	vh = new double[_n];
		
	n = _n;
}

void Vector::Resize(int _n){
	n = _n;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Matrix::Matrix(){
	m  = 0;
	n  = 0;
	l  = 0;
	vh = 0;
}

Matrix::Matrix(int _m, int _n){
	m  = 0;
	n  = 0;
	l  = 0;
	vh = 0;
	Allocate(_m, _n);
}

Matrix::~Matrix(){
	//if(vh) delete[] vh;
}

void Matrix::Delete(){
	if(vh) delete[] vh;
}

Matrix Matrix::SubMatrix(int row, int col, int _m, int _n){
	Matrix y;
	y.vh = vh + row + l*col;
	y.m = _m;
	y.n = _n;
	y.l =  l;
	return y;
}

void Matrix::Allocate(int _m, int _n){
	Delete();
	vh = new double[_m*_n];
	
	m = _m;
	n = _n;
	l = _m;
}

void Matrix::Resize(int _m, int _n){
	m = _m;
	n = _n;
}

ostream& operator<<(ostream& os, Vector& v){
	for(int i = 0; i < v.n; i++)
		os << v(i) << " ";
	return os;
}

ostream& operator<<(ostream& os, Matrix& m){
	for(int i = 0; i < m.m; i++)for(int j = 0; j < m.n; j++)
		os << m(i,j) << " ";
	return os;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void SparseVector::Resize(int _n, int _dim){
	n   = _n;
	dim = _dim;
	//index.resize(n);
	//fill(index.begin(), index.end(), end());
}

const Vector* SparseVector::Item(int i) const{
	//if(index[i] == end())
    const_iterator it = find(i);
    if(it == end())
		return 0;

	//return &(*index[i]).second;
    return &(it->second);
}

Vector SparseVector::SubVector(int i){
    iterator it = find(i);
	//if(index[i] == end()){
    if(it == end()){
		Vector v;
		v.Allocate(dim);
		vec_clear(v);
        insert(make_pair(i, v));
		//index[i] = insert(make_pair(i, v)).first;
		return v;
	}
	//return (*index[i]).second;
    return it->second;
}

void SparseVector::Clear(int i){
    iterator it = find(i);
	//if(index[i] != end()){
    if(it != end()){
        erase(it);
		//erase(index[i]);
		//index[i] = end();
	}
}

void SparseVector::Clear(){
	clear();
	//fill(index.begin(), index.end(), end());
}

///////////////////////////////////////////////////////////////////////////////////////////////////

SparseMatrix::SparseMatrix(){
	//base = 0;
	m   = 0;
	n   = 0;
	dim = 0;
}

void SparseMatrix::Resize(int _m, int _n, int _dim){
	m   = _m;
	n   = _n;
	dim = _dim;
	rows.resize(m);
	//for(Row& r :rows){
	//	r.index.resize(n);
	//	fill(r.index.begin(), r.index.end(), r.end());
	//}
}

const Matrix* SparseMatrix::Item(int i, int j) const{
    Row::const_iterator it = rows[i].find(j);
	//if(rows[i].index[j] == rows[i].end())
    if(it == rows[i].end())
		return 0;

    return &(it->second);
	//return &(*rows[i].index[j]).second;
}

Matrix SparseMatrix::SubMatrix(int i, int j){
    Row::iterator it = rows[i].find(j);
    if(it == rows[i].end()){
	//if(rows[i].index[j] == rows[i].end()){
		Matrix m;
		m.Allocate(dim, dim);
		mat_clear(m);
        rows[i].insert(make_pair(j, m));
		//rows[i].index[j] = rows[i].insert(make_pair(j, m)).first;
		return m;
	}
    return it->second;
	//return (*rows[i].index[j]).second;
}

Matrix SparseMatrix::SubMatrix(int i, int j) const{
    return rows[i].find(j)->second;
	//return (*rows[i].index[j]).second;
}

SparseMatrix SparseMatrix::RowMatrix(int i) const{
	SparseMatrix y;
	y.Resize(1, n, dim);
	for(Row::const_iterator it = rows[i].begin(); it != rows[i].end(); it++){
		int           j  = it->first;
		const Matrix& _m = it->second;
        y.rows[0].insert(make_pair(j, _m));
		//y.rows[0].index[j] = y.rows[0].insert(make_pair(j, _m)).first;
	}
	return y;
}

//SparseMatrix SparseMatrix::ColMatrix(int j) const{
//	SparseMatrix y;
//	y.dim        = dim;
//	y.base       = (SparseMatrix*)this;
//	y.row_offset = 0;
//	y.col_offset = j;
//	
//	for(const_iterator it = begin(); it != end(); it++){
//		if(it->first.second == j){
//			y.insert(make_pair(make_pair(it->first.first, 0), it->second));
//		}
//	}
//	return y;
//}

void SparseMatrix::Clear(){
	for(Row& row : rows){
		row.clear();
		//fill(row.index.begin(), row.index.end(), row.end());
	}
}

void SparseMatrix::RowClear(int i){
	rows[i].clear();
	//fill(rows[i].index.begin(), rows[i].index.end(), rows[i].end());
}

void SparseMatrix::ColClear(int j){
	for(int i = 0; i < m; i++){
		for(Row::iterator it = rows[i].begin(); it != rows[i].end(); ){
			if(it->first == j){
				 it = rows[i].erase(it);
				 //rows[i].index[j] = rows[i].end();
			}
			else it++;
		}
	}
}

void SparseMatrix::PrintSparsity(ostream& os, SparseMatrix::PrintCallback* callback){
	vector<int> ordering(m);
	for(int i = 0; i < m; i++)
		ordering[i] = i;

	if(callback)
		sort(ordering.begin(), ordering.end(), [callback](int lhs, int rhs){ return callback->Order(lhs, rhs); });

	vector<int> ordering_rev(ordering.size());
	for(int i = 0; i < ordering.size(); i++){
		ordering_rev[ordering[i]] = i;
	}
	
	vector<int> col;
	for(int i = 0; i < m; i++){
		int i_ordered = ordering[i];

		col.clear();
		for(Row::iterator it = rows[i_ordered].begin(); it != rows[i_ordered].end(); it++){
			col.push_back(ordering_rev[it->first]);
		}
		sort(col.begin(), col.end());
		
		int j = 0;
		for(vector<int>::iterator it = col.begin(); it != col.end(); it++){
			while(j < *it){
				os << (callback ? callback->Label(ordering[i], ordering[j], false) : '0') << ' ';
				j++;
			}
			os << (callback ? callback->Label(ordering[i], ordering[j], true) : '1') << ' ';
			j++;
		}
		while(j < n){
			os << (callback ? callback->Label(ordering[i], ordering[j], false) : '0') << ' ';
			j++;
		}
		os << '\n';
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

LinearSolverCustom::LinearSolverCustom(){
	mindeg = true;
}

void LinearSolverCustom::Init(SparseMatrix& A){
	if(A.m == 0)
		return;

	Aii         .resize(A.m);
	Aii_inv     .resize(A.m);
	Arow        .resize(A.m);
	bi          .resize(A.m);
	Aii_inv_bi  .resize(A.m);
	Aii_inv_Arow.resize(A.m);

	for(int i = 0; i < A.m; i++){
		Aii         [i].Allocate(A.dim, A.dim);
		Aii_inv     [i].Allocate(A.dim, A.dim);
		Arow        [i].Resize(1, A.n, A.dim);
		bi          [i].Allocate(A.dim);
		Aii_inv_bi  [i].Allocate(A.dim);
		Aii_inv_Arow[i].Resize(1, A.n, A.dim);
	}
}

void LinearSolverCustom::Finish(){

}

void LinearSolverCustom::Solve(SparseMatrix& A, SparseVector& b, SparseVector& x){
	if(A.m == 0)
		return;

	for(int i = 0; i < A.m; i++){
		Arow        [i].Clear();
		Aii_inv_Arow[i].Clear();
	}

	queue.clear();
	order.clear();

	for(int i = 0; i < A.m; i++){
		queue.insert(i);
	}

	while(!queue.empty()){
		// index with minimum degree
		//int i = A.MinimumDegreeIndex();
		set<int>::iterator it;
		if(mindeg){
			it = min_element(queue.begin(), queue.end(), [&A](int i0, int i1){ return A.rows[i0].size() < A.rows[i1].size(); });
		}
		else{
			it = queue.begin();
		}
		int i = *it;
		order.push_back(i);
		queue.erase(it);

		mat_copy(A.RowMatrix(i), Arow[i]);
		mat_copy(Arow[i].SubMatrix(0, i), Aii[i]);
		vec_copy(b.SubVector(i), bi[i]);

		Arow[i].ColClear(i);
		A.ColClear(i);
		A.RowClear(i);
		b.Clear(i);

		mat_inv_pd(Aii[i], Aii_inv[i]);

		symmat_vec_mul(Aii_inv[i], bi[i], Aii_inv_bi[i], 1.0, 0.0);

		mattr_vec_mul(Arow[i], Aii_inv_bi[i], b, -1.0, 1.0);

		symmat_mat_mul(Aii_inv[i], Arow[i]        , Aii_inv_Arow[i],  1.0, 0.0);
		mattr_mat_mul (Arow[i]   , Aii_inv_Arow[i], A              , -1.0, 1.0);
	}

	for(vector<int>::reverse_iterator it = order.rbegin(); it != order.rend(); it++){
		int i = *it;
		mat_vec_mul   (Arow[i]   , x    , bi[i]         ,  1.0, 1.0);
		
		Vector xi = x.SubVector(i);
		symmat_vec_mul(Aii_inv[i], bi[i], xi, -1.0, 0.0);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_CHOLMOD
# include <cholmod.h>
#endif

struct CholmodWork{
#ifdef USE_CHOLMOD
	cholmod_common   cm;
	cholmod_triplet* A_trip;
	cholmod_sparse*  A_sparse;
	cholmod_factor*  L;
	cholmod_dense*   b_dense;
	cholmod_dense*   x_dense;
#endif
};

LinearSolverCholmod::LinearSolverCholmod(){

}

void LinearSolverCholmod::Init(SparseMatrix& A){
	if(A.m == 0)
		return;

#ifdef USE_CHOLMOD
	nrow = A.m;
	ncol = A.n;
	nblock = 0;

	for(auto& row : A.rows)
		nblock += (int)row.size();

	nzmax = nblock*(A.dim*A.dim);

	CholmodWork* w = new CholmodWork();
	cholmod_start(&w->cm);

	w->cm.nmethods           = 1;
	w->cm.method[0].ordering = CHOLMOD_AMD;
	w->cm.supernodal         = CHOLMOD_AUTO; //CHOLMOD_SIMPLICIAL;
	
	w->A_trip  = cholmod_allocate_triplet(A.dim*nrow, A.dim*ncol, nzmax, 1, CHOLMOD_REAL, &w->cm);
	w->b_dense = cholmod_zeros(A.dim*nrow, 1, CHOLMOD_REAL, &w->cm);

	work = w;
#endif
}

void LinearSolverCholmod::Finish(){
#ifdef USE_CHOLMOD
	CholmodWork* w = (CholmodWork*)work;

	cholmod_free_dense  (&w->b_dense , &w->cm);
	cholmod_free_triplet(&w->A_trip  , &w->cm);

	cholmod_finish(&w->cm);

	delete (CholmodWork*)work;
#endif
}

void LinearSolverCholmod::Solve(SparseMatrix& A, SparseVector& b, SparseVector& x){
	if(A.m == 0)
		return;

#ifdef USE_CHOLMOD	
	CholmodWork* w = (CholmodWork*)work;

	int idx = 0;
	for(int ix0 = 0; ix0 < A.rows.size(); ix0++){
		for(auto it = A.rows[ix0].begin(); it != A.rows[ix0].end(); it++){
			int     ix1 = it->first;
			Matrix& _m  = it->second;

			if(ix0 > ix1)
				continue;

			for(int i = 0; i < A.dim; i++){
				int j0 = (ix0 == ix1 ? i : 0);
				for(int j = j0; j < A.dim; j++){
					((int   *)w->A_trip->i)[idx] = A.dim*ix0 + i;
					((int   *)w->A_trip->j)[idx] = A.dim*ix1 + j;
					((double*)w->A_trip->x)[idx] = _m(i,j);
					idx++;
				}
			}
		}
	}
	w->A_trip->nnz = idx;

	for(auto it = b.begin(); it != b.end(); it++){
		int     ix = it->first;
		Vector& _v = it->second;

		for(int i = 0; i < A.dim; i++){
			((double*)w->b_dense->x)[A.dim*ix + i] = _v(i);
		}
	}

	w->A_sparse = cholmod_triplet_to_sparse(w->A_trip, 0, &w->cm);

	w->L = cholmod_analyze (w->A_sparse, &w->cm);
	cholmod_factorize (w->A_sparse, w->L, &w->cm) ;
	w->x_dense = cholmod_solve (CHOLMOD_A, w->L, w->b_dense, &w->cm);

	// store result
	for(int ix = 0; ix < ncol; ix++){
		Vector& _v = x.SubVector(ix);
		for(int j = 0; j < A.dim; j++)
			_v(j) = -1.0 * ((double*)w->x_dense->x)[A.dim*ix + j];  //< flip sign here
	}

	cholmod_free_sparse (&w->A_sparse, &w->cm);
	cholmod_free_factor (&w->L       , &w->cm);
	cholmod_free_dense  (&w->x_dense , &w->cm);
#endif
}

}
