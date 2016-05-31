#ifndef ap_utils_math_tpp__
#define ap_utils_math_tpp__

//#include <algorithm>
//#include <cassert>
#include <numeric>

//------------------------------------------------------------------------------------------
template <typename t> void
axblu(t* const a, double* x, t* const b, std::size_t const m, std::size_t const p){
  // solve for X in AX=B by LU decomposition method. 
  // m = numRows A
  // m = numCols A
  // p = numCols X

  std::vector<double> aa(a, a+m*m);
  gsl_matrix_view lu = gsl_matrix_view_array(aa.data(), m, m);

  // result will be stored in this
  gsl_matrix_view gx = gsl_matrix_view_array(x, m, p);
  
  int s;
  gsl_permutation* perm = gsl_permutation_alloc (m);

  assert(GSL_SUCCESS == gsl_linalg_LU_decomp(&lu.matrix, perm, &s));
  
  gsl_matrix_view gb = gsl_matrix_view_array(b, m, p);

  gsl_vector_view btemp;
  gsl_vector_view xtemp;
  for(std::size_t col = 0; col < p; ++col){
    btemp = gsl_matrix_column(&gb.matrix, col);
    xtemp = gsl_matrix_column(&gx.matrix, col);
    gsl_linalg_LU_solve(&lu.matrix, perm, &btemp.vector, &xtemp.vector);
  }

  gsl_permutation_free(perm); 

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t> void
axbsvd(double * x, const t * a, const t * b, size_t const m, size_t const n, size_t const p){
  // solve for X in AX=B by SVD decomposition method. 
  // m == numRows A == numRows B
  // n == numCols A == numRows X
  // p == numCols X == numCols B

  gsl_matrix * U = gsl_matrix_alloc(m, n);
  gsl_matrix * V = gsl_matrix_alloc(n, n);
  gsl_vector * s = gsl_vector_alloc(n);
  gsl_vector * w = gsl_vector_alloc(n);

  std::copy(a, a+m*n, U->data);
  assert(GSL_SUCCESS == gsl_linalg_SV_decomp(U, V, s, w));
  
  for(size_t i = 0; i < p; ++i){
    gsl_vector_view xcol = gsl_vector_view_array_with_stride(x+i, p, n);
    gsl_vector_const_view bcol = gsl_vector_const_view_array_with_stride(b+i, p, m);
    gsl_linalg_SV_solve(U, V, s, &bcol.vector, &xcol.vector);
  }
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(s);
  gsl_vector_free(w);

  return;
}
//------------------------------


///////////////////////////////////////////////////////////////////////////////
template<typename t> t*
eye(std::size_t const n){
  // identity matrix

  //http://stackoverflow.com/questions/808464/c-new-call-that-behaves-like-calloc
  t* tmp = new t[n*n](); // value initialization
  for(std::size_t i = 0; i < n; ++i) tmp[i*(n+1)] = 1;
  return tmp;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t> std::vector <t>
subtract(const std::vector <t> & a, const std::vector <t> & b){
  size_t const size = a.size();
  assert(size==b.size());
  std::vector<t> result(size);
  std::transform(a.cbegin(), a.cend(), b.cbegin(), result.begin(), std::minus<t>());
  return result;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <class t> std::vector <t>
add(const std::vector <t> & a, const std::vector <t> & b){
  size_t const size = a.size();
  assert(size==b.size());
  std::vector<t> result(size);
  std::transform(a.cbegin(), a.cend(), b.cbegin(), result.begin(), std::plus<t>());
  return result;
}
//------------------------------


///////////////////////////////////////////////////////////////////////////////
template <typename t> t sum(const t * base, size_t const size){
  assert(base);
  t ans = 0;
  for(size_t i = 0; i < size; ++i) ans += *(base+i);
  return ans;
}
//------------------------------

//------------------------------------------------------------------------------------------
template <typename t> t
prod(const t * base, size_t const size){
  assert(base);
  t ans = 1;
  for(size_t i = 0; i < size; ++i) ans *= *(base+i);
  return ans;
}
//------------------------------

// template <typename input_iterator>
// typename input_iterator::value_type prod(input_iterator begin, input_iterator end) {
//   return std::accumulate(begin, end, 1, std::multiplies<typename input_iterator::value_type>());
// }


template <typename I>
typename std::iterator_traits<I>::value_type prod(I begin, I end) {
  return std::accumulate(begin, end, typename std::iterator_traits<I>::value_type(1),
                         std::multiplies<typename std::iterator_traits<I>::value_type>());
}


//------------------------------------------------------------------------------------------
template <typename t> double
one_norm(const std::vector<t>& array){
  return cblas_dasum(array.size(), array.data(), 1);
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t> double
twoNorm(const std::vector<t>& array){
  return cblas_dnrm2(array.size(), array.data(), 1);
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t> double
twoNorm(const t * base1, const t * base2, size_t size){
  // twoNorm of (vec0-vec1)
  assert(base1);
  assert(base2);
  std::vector<double> diff(size);
  for(size_t n = 0; n < size; ++n){
    diff[n] = (base1[n] - base2[n]);
  }
  return cblas_dnrm2(size, diff.data(), 1);
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t> double
twoNorm(const std::vector<t>& vec0, const std::vector<t>& vec1){
  // twoNorm of (vec0-vec1)
  
  size_t const size = vec0.size();
  assert(size == vec1.size());
  std::vector<double> diff(size);
  for(size_t n = 0; n < size; ++n){
    diff[n] = (vec0[n] - vec1[n]);
  }
  return cblas_dnrm2(size, diff.data(), 1);
}
//------------------------------

//------------------------------------------------------------------------------------------
template <typename t> size_t
rank(const t * base, size_t const nrows, size_t const ncols, double const tol){
  auto u = gsl_matrix_alloc(nrows, ncols);
  auto v = gsl_matrix_alloc(ncols, ncols);
  auto s = gsl_vector_alloc(ncols);

  auto a = gsl_matrix_const_view_array(base, nrows, ncols);
  gsl_matrix_memcpy(u, &a.matrix);

  assert(GSL_SUCCESS == gsl_linalg_SV_decomp_jacobi(u, v, s));
  size_t rank = 0;
  for(size_t i = 0; i < ncols; ++i){
    if(tol < gsl_vector_get(s, i)){
      ++rank;
    }
  }

  gsl_matrix_free(u);
  gsl_matrix_free(v);
  gsl_vector_free(s);
  
  return rank;
}
//------------------------------




#endif 
