#ifndef __imoga_stats_tpp__
#define __imoga_stats_tpp__

#include <gsl/gsl_statistics.h>

//------------------------------------------------------------------------------------------
template <typename t> double
mean(t* const data, std::size_t const stride, std::size_t const size)
{
  double m = 0.0;
  for(std::size_t i = 0; i < size; ++i){m += (data[i*stride] - m)/(i+1);}
  
  return m;
}
//------------------------------


// //------------------------------------------------------------------------------------------
// template <typename t> void
// mean(double* means, t* const matrix, std::size_t const nrows, std::size_t const ncols)
// {
//   std::fill(means, means+ncols, 0.0);
  
//   for(std::size_t i = 0; i < nrows; ++i)
//     for(std::size_t j = 0; j < ncols; ++j)
//       {
// 	means[j] += matrix[i*ncols+j];
//       }

//   for(std::size_t i = 0; i < ncols; ++i) means[i] /= nrows;
// }
// //------------------------------


//------------------------------------------------------------------------------------------
template <typename t> void
mean(double* means, t* const matrix, std::size_t const nrows, std::size_t const ncols)
{
  //std::cout << "\n"__FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
  for(std::size_t j = 0; j < ncols; ++j)
    {
      means[j] = mean(matrix+j, ncols, nrows);
    }
}
//------------------------------


//------------------------------------------------------------------------------------------
template<typename t> void
var(double* vars, t* const matrix, std::size_t const nrows, std::size_t const ncols)
{
  //https://secure.wikimedia.org/wikipedia/en/wiki/Algorithms_for_calculating_variance
  
  for(std::size_t j = 0; j < ncols; ++j)
    {
      vars[j] = gsl_stats_variance(matrix+j, ncols, nrows);
    }
}
//------------------------------


//------------------------------------------------------------------------------------------
template<typename t> void
stdDev(double* s, t* const matrix, std::size_t const nrows, std::size_t const ncols)
{
  var(s, matrix, nrows, ncols);
  std::for_each(s, s+ncols, [](double&i){ i = sqrt(i); });
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t> void
studentize(double* out, t* const in, std::size_t const nrows, std::size_t const ncols)
{
  double* means = new double [ncols];
  mean(means, in, nrows, ncols);

  double* sd = new double [ncols];
  stdDev(sd, in, nrows, ncols);

  for(std::size_t i = 0; i < nrows; ++i)
    for(std::size_t j = 0; j < ncols; ++j)
      out[i*ncols+j] = (in[i*ncols+j] - means[j])/sd[j];

  delete [] means;
  delete [] sd;
}
//------------------------------


//------------------------------------------------------------------------------------------
template<typename t> void
center(double* out, t* const in, std::size_t const nrows, std::size_t const ncols, t* const means)
{
  for(std::size_t i = 0; i < nrows; ++i)
    for(std::size_t j = 0; j < ncols; ++j)
      out[i*ncols+j] = (in[i*ncols+j] - means[j]);
}
//------------------------------


//------------------------------------------------------------------------------------------
template<typename t> void
center(t* matrix, std::size_t const nrows, std::size_t const ncols)
{
  double* means = new double [ncols];
  
  mean(means, matrix, nrows, ncols);

  for(std::size_t i = 0; i < nrows; ++i)
    for(std::size_t j = 0; j < ncols; ++j)
      matrix[i*ncols+j] -= means[j];

  delete [] means;
}
//------------------------------


#endif
