#ifndef __AP_UTILS_GSL_H__
#define __AP_UTILS_GSL_H__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include <cstddef>
#include <cassert>
#include <string>
#include <iosfwd>

enum svd_method {normal, mod, jacobi};

#define AP_GSL_SHOW(x) std::cout << "\n" #x "\n"; ap_gsl_show((x))
#define AP_GSL_WRITE(x) ap_gsl_write(x, #x".out")

inline std::size_t ap_gsl_matrix_numel(gsl_matrix * const m){
  assert(m);
  return m->size1*m->size2;
}

void ap_gsl_set_random(gsl_matrix * m, const gsl_rng * r);

void ap_gsl_show(gsl_matrix const * m);
void ap_gsl_show(gsl_vector const * m);

void ap_gsl_write(gsl_matrix const * m, const std::string & fileName);
void ap_gsl_write(gsl_vector const * v, const std::string & fileName);

double ap_gsl_dasum(gsl_matrix const * m1, const gsl_matrix * m2);
double ap_gsl_dasum(gsl_matrix const * m);

double ap_gsl_sum(gsl_vector const * v);
double ap_gsl_prod(gsl_vector const * v);


int ap_svd(gsl_matrix * u, gsl_vector * s, gsl_matrix * v,
           gsl_matrix const * a);
int ap_svd_jacobi(gsl_matrix * u, gsl_vector * s, gsl_matrix * v,
                  gsl_matrix const * a);

int ap_gsl_linalg_SV_solve(gsl_matrix * x, gsl_matrix const * u,
                           gsl_vector const * s, gsl_matrix const * v,
                           gsl_matrix const * b);

inline void ap_gsl_free(gsl_matrix * m){
  gsl_matrix_free(m);
}

inline void ap_gsl_free(gsl_vector * m){
  gsl_vector_free(m);
}

gsl_matrix * ap_gsl_clone(gsl_matrix const * m);
gsl_vector * ap_gsl_clone(gsl_vector const * v);

void ap_singular_values(gsl_matrix const * m);

#endif /* AP_UTILS_GSL_H__ */
