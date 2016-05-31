#include "ap_utils_gsl.h"
#include "ap_utils.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <limits>
#include <vector>
#include <iostream>
#include <fstream>

gsl_matrix* ap_gsl_clone(gsl_matrix const * m){
  gsl_matrix * clone = gsl_matrix_alloc(m->size1, m->size2);
  gsl_matrix_memcpy(clone, m);
  return clone;
}

gsl_vector* ap_gsl_clone(gsl_vector const * v){
  gsl_vector * clone = gsl_vector_alloc(v->size);
  gsl_vector_memcpy(clone, v);
  return clone;
}

void ap_gsl_set_random(gsl_matrix * m, const gsl_rng * r){
  for(size_t i = 0; i < m->size1; ++i){
    for(size_t j = 0; j < m->size2; ++j){
      gsl_matrix_set(m, i, j, gsl_rng_uniform(r));
    }
  }
}


void ap_gsl_show(gsl_matrix const * m){
  std::cout << "# " << m->size1 << " x " << m->size2 << " :: " << __FILE__ << " --> "
	    << __func__ << "() --> " << __LINE__ << "\n";
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  for(size_t i = 0; i < m->size1; ++i){
    for(size_t j = 0; j < m->size2; ++j){
      std::cout << gsl_matrix_get(m, i, j) << "\t";
    }
    std::cout << "\n";
  }
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
}


void ap_gsl_show(gsl_vector const * v){
  std::cout << "# " << v->size << " elements :: " << __FILE__ << " --> "
	    << __func__ << "() --> " << __LINE__ << "\n";
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  for(size_t i = 0; i < v->size; ++i){
    std::cout << gsl_vector_get(v, i) << "\t";
  }
  std::cout << "\n";
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
}


void ap_gsl_copy_data(double * v, const gsl_matrix * m){
  for(size_t i = 0, n = 0; i < m->size1; ++i){
    for(size_t j = 0; j < m->size2; ++j, ++n){
      v[n] = gsl_matrix_get(m, i, j);
    }
  }
}
  

double ap_gsl_dasum(const gsl_matrix * m){
  std::vector<double> v(m->size1*m->size2);
  ap_gsl_copy_data(v.data(), m);
  return cblas_dnrm2(v.size(), v.data(), 1);
}
  

double ap_gsl_dasum(const gsl_matrix * m1, const gsl_matrix * m2){
  gsl_matrix * diff = gsl_matrix_alloc(m1->size1, m1->size2);
  gsl_matrix_memcpy(diff, m1);
  gsl_matrix_sub(diff, m2);
  double dasum = ap_gsl_dasum(diff);
  gsl_matrix_free(diff);
  return dasum;
}


double ap_gsl_sum(const gsl_vector *v){
  double sum{0.0};
  for(size_t i = 0; i < v->size; ++i) sum += gsl_vector_get(v, i);
  return sum;
}


double ap_gsl_prod(const gsl_vector *v){
  double prod{1.0};
  for(size_t i = 0; i < v->size; ++i) prod *= gsl_vector_get(v, i);
  return prod;
}


void ap_gsl_write(const gsl_matrix * m, const std::string & fileName){
  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good()){
#ifdef AP_HEADERS
    outFile << "#"<< __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
    outFile << "# " << m->size1 << " x " << m->size2 << "\n";
#endif
      
    outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    outFile.precision(std::numeric_limits<double>::digits10);

    for(size_t row = 0; row < m->size1; ++row){
#ifdef AP_ROW_NUMBERS
      outFile << row << "\t";
#endif
      for(size_t col = 0; col < m->size2; ++col){
        outFile << gsl_matrix_get(m, row, col) << "\t"; 
      }
      outFile << "\n";
    }
    outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint); 
  } else {
    std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
    std::cerr << alert("ERROR: ") << "could not open file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }
}


void ap_gsl_write(const gsl_vector * v, const std::string & fileName){
  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good()){
#ifdef AP_HEADERS
    outFile << "#"<< __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
    outFile << "# " << v->size << "elements\n";
#endif
      
    outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    outFile.precision(std::numeric_limits<double>::digits10);

      for(size_t col = 0; col < v->size; ++col){
        outFile << gsl_vector_get(v, col) << "\t";
      }
      outFile << "\n";

    outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint); 
  } else {
    std::cout << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
    std::cerr << alert("ERROR: ") << "could not open file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  return;
}


//------------------------------------------------------------------------------------------
int svd_solve(gsl_matrix * x, gsl_matrix const * a, gsl_matrix const * b,
              double const tol){
  // solve for X in AX=B by SVD decomposition method. 
  // m == numRows A == numRows B
  // n == numCols A == numRows X
  // p == numCols X == numCols B

  std::cout << "\n" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
  
  const size_t m = a->size1;
  const size_t n = a->size2;
  const size_t p = b->size2;

  gsl_matrix * U = gsl_matrix_alloc(m, n);
  gsl_matrix * V = gsl_matrix_alloc(n, n);
  gsl_vector * s = gsl_vector_alloc(n);
  gsl_vector * w = gsl_vector_alloc(n);

  gsl_matrix_memcpy(U, a);
  
  assert(GSL_SUCCESS == gsl_linalg_SV_decomp(U, V, s, w));

  {
    //double norm = gsl_vector_get(s, 0); // largest singular value
    //SHOW(std::max(m,n)*norm*std::numeric_limits<double>::epsilon());
    //SHOW(gsl_vector_get(s, 0)/gsl_vector_get(s, n-1));
  }

  //AP_GSL_SHOW(s);

  for(size_t i = 0; i < n; ++i){
    if(tol > gsl_vector_get(s, i)){
      gsl_vector_set(s, i, 0.0);
    }
  }

  for(size_t i = 0; i < p; ++i){
    gsl_vector_view xcol = gsl_matrix_column(x,i);
    gsl_vector_const_view bcol = gsl_matrix_const_column(b, i); 
    gsl_linalg_SV_solve(U, V, s, &bcol.vector, &xcol.vector);
  }
  
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(s);
  gsl_vector_free(w);

  return GSL_SUCCESS;
}


///////////////////////////////////////////////////////////////////////////////
int ap_svd(gsl_matrix * u, gsl_vector * s, gsl_matrix * v,
           gsl_matrix const * a){
  const size_t n = a->size2;
  gsl_vector * w = gsl_vector_alloc(n);
  gsl_matrix_memcpy(u, a);
  assert(GSL_SUCCESS == gsl_linalg_SV_decomp(u, v, s, w));
  gsl_vector_free(w);
  return GSL_SUCCESS;
}



///////////////////////////////////////////////////////////////////////////////
int ap_svd_jacobi(gsl_matrix * u, gsl_vector * s, gsl_matrix * v,
                  gsl_matrix const * a){
  gsl_matrix_memcpy(u, a);
  assert(GSL_SUCCESS == gsl_linalg_SV_decomp_jacobi(u, v, s));
  return GSL_SUCCESS;
}


///////////////////////////////////////////////////////////////////////////////
int ap_gsl_linalg_SV_solve(gsl_matrix * x, gsl_matrix const * u,
                           gsl_vector const * s, gsl_matrix const * v,
                           gsl_matrix const * b){
  const size_t p = x->size2;
  for(size_t i = 0; i < p; ++i){
    gsl_vector_view xcol = gsl_matrix_column(x,i);
    gsl_vector_const_view bcol = gsl_matrix_const_column(b, i); 
    gsl_linalg_SV_solve(u, v, s, &bcol.vector, &xcol.vector);
  }
  return GSL_SUCCESS;
}


///////////////////////////////////////////////////////////////////////////////
void ap_singular_values(gsl_matrix const * m){
  auto u = ap_gsl_clone(m);
  auto s = gsl_vector_alloc(m->size2);
  auto w = gsl_vector_alloc(m->size2);
  auto v = gsl_matrix_alloc(m->size2, m->size2);
  gsl_linalg_SV_decomp(u, v, s, w);
  ap_gsl_show(s);
  ap_gsl_free(u);
  ap_gsl_free(s);
  ap_gsl_free(w);
  ap_gsl_free(v);
  return;
}
