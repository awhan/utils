#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath> 
#include <vector>
#include <cfloat> // for DBL_EPSILON
#include <sys/stat.h>
#include <cassert>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h> // for gsl_ran_shuffle()
#include <gsl/gsl_sort.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h> // for gsl_fcmp


#include "ap_utils.h"


///////////////////////////////////////////////////////////////////////////////
size_t rank(gsl_matrix * const m){
  size_t const nrows = m->size1;
  size_t const ncols = m->size2;
  auto u = gsl_matrix_alloc(nrows, ncols);
  auto v = gsl_matrix_alloc(ncols, ncols);
  auto s = gsl_vector_alloc(ncols);
  auto w = gsl_vector_alloc(ncols);
  
  gsl_matrix_memcpy(u, m);
  
  assert(GSL_SUCCESS == gsl_linalg_SV_decomp(u, v, s, w));

  double const tol = std::max(nrows, ncols) * gsl_vector_get(s,0) *
    std::numeric_limits<double>::epsilon();
  size_t rank = 0;
  for(size_t i = 0; i < ncols; ++i){
    if(tol < gsl_vector_get(s, i)){
      ++rank;
    }
  }
  
  gsl_matrix_free(u);
  gsl_matrix_free(v);
  gsl_vector_free(s);
  gsl_vector_free(w);
  return rank;
}
//------------------------------


//------------------------------------------------------------------------------------------
void center(gsl_matrix* matrix, gsl_vector* const means)
{
  display(matrix->data, matrix->size1, matrix->size2);
  gsl_vector_view v;
  for(size_t row = 0; row < matrix->size1; ++row)
    {
      v = gsl_matrix_row(matrix, row);
      gsl_vector_sub(&v.vector, means);
    }
      
}
//------------------------------

//------------------------------------------------------------------------------------------
void center(std::vector<std::vector<double> >& matrix, const int flag)
{
  if(1 == flag)
    {
      std::vector<double> m(matrix[0].size());
      mean(m, matrix, 1);
      for(auto &i: matrix)
	{
	  std::transform(i.begin(), i.end(), m.begin(), i.begin(),std::minus<double>());
	}
    }
  else if(2 == flag)
    {
      std::vector<double> m(matrix[0].size());
      mean(m, matrix, 1);
      for(auto &i: matrix)
	{
	  std::transform(i.begin(), i.end(), m.begin(), i.begin(),std::minus<double>());
	}
    }
  else
    {
      std::cerr << "\n" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      std::cout << redColor << "ERROR" << noColor << std::endl;
    }
}
//------------------------------

//------------------------------------------------------------------------------------------
void mean(std::vector<double>& means, const std::vector<std::vector<double> >& matrix,
	  const int flag)
{
  for(auto &i: means){i = 0.0;}
  
  if(1 == flag)//column means
    {
      for(auto &i: matrix)
	{
	  std::transform(means.begin(), means.end(), i.begin(), means.begin(),
			 std::plus<double>());
	}
      const unsigned numRows = matrix.size();
      for(auto &i: means){ i /= numRows;}
    }
  else if(2 == flag)//row means
    {
      std::transform(matrix.begin(), matrix.end(), means.begin(),
		     [](const std::vector<double>& i)
      		     {
      		       return std::accumulate(i.begin(), i.end(), 0.0)/i.size();
      		     });
    }
  else
    {
      std::cerr << "\n" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      std::cout << redColor << "ERROR" << noColor << std::endl;
    }
}
//------------------------------

//------------------------------------------------------------------------------------------
void mean(gsl_vector* means, const gsl_matrix* matrix)
{
  for(size_t col = 0; col < matrix->size2; ++col)
    {
      gsl_vector_const_view v = gsl_matrix_const_column(matrix, col);
      gsl_vector_set(means, col, gsl_stats_mean(v.vector.data, v.vector.stride, v.vector.size));
    }
}
//------------------------------


//------------------------------------------------------------------------------------------
void distanceMatrix(std::vector<std::vector<double> >& dm,
		    const std::vector<std::vector<double> >& data,
		    const std::vector<std::vector<double> >& centers)
{
  const unsigned numData = data.size();
  const unsigned numCenters = centers.size();
  for(unsigned i = 0; i < numData; ++i)
    {
      for(unsigned j = 0; j < numCenters; ++j)
	{
	  dm[i][j] = twoNorm<double>(data[i], centers[j]);
	}
    }
}
//------------------------------

//------------------------------------------------------------------------------------------
void
distanceMatrix(double* out, double* const data, size_t const numData, double* const centers,
		size_t const numCenters, size_t const ndims)
{
  //ASSUMPTION:: data is a matrix of numData rows and ndims columns,
  //centers is a matrix of numCenters rows and ndims columns

  double* temp = new double [ndims];
  
  for(size_t i = 0; i < numData; ++i)
    {
      for(size_t j = 0; j < numCenters; ++j)
	{
	  std::copy(data+i*ndims, data+i*ndims+ndims, temp); 
	  cblas_daxpy(ndims, -1.0, centers+j*ndims, 1, temp, 1);
	  out[i*numCenters+j] = cblas_dnrm2(ndims, temp, 1);
	}
    }

  delete [] temp;
}
//------------------------------

// //------------------------------------------------------------------------------------------
// bool nonDominated(const std::vector<double>& v, const std::vector<std::vector<double> >& a) 
// {
//   // is f non-dominated?
//   for(auto itr = a.begin(); itr != a.end(); ++itr)
//     {
//       if(dominanceCompare(v, *itr)) return false;
//     }
//   return true;
// }
// //------------------------------

//------------------------------------------------------------------------------------------
//Thu Sep 16 18:55:35 2010
void pinv(gsl_matrix* apinv, const gsl_matrix* a)
{
  // numRows
  const unsigned  m = a->size1;

  // numCols 
  const unsigned  n = a->size2;
  
  assert(m == apinv->size2);
  assert(n == apinv->size1);
  
  gsl_matrix* u = gsl_matrix_alloc(m, n);
  gsl_matrix_memcpy(u, a);

  // for singular values
  gsl_vector* s = gsl_vector_alloc(n);

  // right singular vectors
  gsl_matrix* v = gsl_matrix_alloc(n, n);

  gsl_vector* w = gsl_vector_alloc(n);

  // std::cout << "before svd" << std::endl;
  // std::cout << "u" << std::endl;
  // display<double>(u->data, m, n);

  gsl_linalg_SV_decomp(u, v, s, w);

  // std::cout << "after svd" << std::endl;
  // std::cout << "u" << std::endl;
  // display<double>(u->data, m, n);

  // std::cout << "v" << std::endl;
  // display<double>(v->data, n, n);
  
  // std::cout << "singular values" << std::endl;
  // display<double>(s->data, n);

  // form the pseudo-inverse

  // form product of V and \Sigma^
  for(unsigned col = 0; col < n; ++col)
    {
      gsl_vector_view temp = gsl_matrix_column(v, col);
      gsl_vector_scale(&temp.vector, 1.0/s->data[col]);
    }

  // std::cout << "v matrix after scaling is" << std::endl;
  // display<double>(v->data, n, n);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, v, u, 0.0, apinv);

  // std::cout << "pseudo-inverse of A is" << std::endl;
  // display<double>(apinv->data, n, m);

  gsl_matrix_free(u);
  gsl_matrix_free(v);
  gsl_vector_free(s);
  gsl_vector_free(w);
  
  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
void axbpinv(const std::vector<std::vector<double> >& a,
	     std::vector<std::vector<double> >& x,
	     const std::vector<std::vector<double> >& b,
	     const unsigned & m, const unsigned & n, const unsigned & p)
{
  //Thu Sep 16 18:55:23 2010
  // TODO return the residual as well
  // m = numRows A
  // n = numCols A
  // p = numCols X
  
  gsl_matrix* apinv = gsl_matrix_calloc(n, m);
  gsl_matrix* ga = gsl_matrix_alloc(m, n);
  stl2gsl(a, ga, m, n);

  if(n>m) 
    {// fat matrix
      gsl_matrix* aTrans = gsl_matrix_calloc(n, m);
      gsl_matrix_transpose_memcpy(aTrans, ga);
      
      // pinv of A transpose
      gsl_matrix* aTransPinv = gsl_matrix_alloc(m, n);
      pinv(aTransPinv, aTrans);

      // apinv = transpose(aTransPinv) 
      gsl_matrix_transpose_memcpy(apinv, aTransPinv);

      gsl_matrix_free(aTrans);
      gsl_matrix_free(aTransPinv);	  
    }
  else
    {// thin or square matrix
      pinv(apinv, ga);
    }

  gsl_matrix* gx = gsl_matrix_alloc(n, p);
  gsl_matrix* gb = gsl_matrix_alloc(m, p);
  stl2gsl(b, gb, m, p);
  
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, apinv, gb, 0.0, gx);

  gsl2stl(gx, x, n, p);


  gsl_matrix_free(ga);
  gsl_matrix_free(apinv);
  gsl_matrix_free(gx);
  gsl_matrix_free(gb);

  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
void axblu(const std::vector<std::vector<double> >& a, std::vector<std::vector<double> >& x,
	   const std::vector<std::vector<double> >& b, size_t const m, size_t const p){
  // solve for X in AX=B by LU decomposition method. 
  // m = numRows A
  // m = numCols A
  // p = numCols X

  // copy is made as lu is decomposed in place
  gsl_matrix* lu = gsl_matrix_alloc(m, m);
  stl2gsl(a, lu, m, m);

  // result will be stored in this
  gsl_matrix* gx = gsl_matrix_alloc(m, p);
  
  int s;
  gsl_permutation* perm = gsl_permutation_alloc (m);

  assert(GSL_SUCCESS == gsl_linalg_LU_decomp(lu, perm, &s));
  
  gsl_matrix* gb = gsl_matrix_alloc(m, p);
  stl2gsl(b, gb, m, p);

  gsl_vector_view btemp;
  gsl_vector_view xtemp;
  for(unsigned col = 0; col < p; ++col)
    {
      btemp = gsl_matrix_column(gb, col);
      xtemp = gsl_matrix_column(gx, col);
      gsl_linalg_LU_solve(lu, perm, &btemp.vector, &xtemp.vector);
    }

  gsl2stl(gx, x, m, p);
  
  gsl_matrix_free(lu);
  gsl_matrix_free(gx);
  gsl_matrix_free(gb);
  gsl_permutation_free(perm); 

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
void stl2gsl(const std::vector<std::vector<double> >& stl, gsl_matrix* gsl,
	     const unsigned & numRows, const unsigned & numCols)
{
  // convert STL in to GSL
  
  for(unsigned row = 0; row < numRows; ++row)
    {
      for(unsigned col =0; col < numCols; ++col)
	{
	  gsl->data[row*numCols+col] = stl[row][col];
	}
    }
  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
void stl2gsl(const std::vector<double>& stl, gsl_vector* gsl)
{
  // convert STL in to GSL
  const unsigned numElements = stl.size();

  gsl = gsl_vector_alloc(numElements);

  std::copy(stl.begin(), stl.end(), gsl->data);
  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
void gsl2stl(const gsl_matrix* gsl, std::vector<std::vector<double> >& stl,
	     const unsigned & numRows, const unsigned & numCols)
{
  // convert GSL in to STL 
  for(unsigned row = 0; row < numRows; ++row)
    {
      for(unsigned col =0; col < numCols; ++col)
	{
	  stl[row][col] = gsl->data[row*numCols+col];
	}
    }
  return;
}
//------------------------------  






//------------------------------------------------------------------------------------------
void makeFile(std::ofstream & ofs, const std::vector<std::string>& info)
{
  const std::string fileName = info[0];
  const std::string fileCreatedIn = info[1];
  const std::string functionCreatedIn = info[2];
  const std::string lineNumber = info[3];

  if (fileExists(fileName))
    {
      std::cerr << __FILE__ << " --> " << __FUNCTION__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR:: ") << fileName 
		<< " already exists; won't overwrite ... " << std::endl;
      exit(EXIT_FAILURE);
    }

  ofs.open(fileName.c_str(), std::ios::out);

  if(ofs.good())
    {
      ofs << "# written by :" << __FILE__ << " --> " << __FUNCTION__ 
	      << "() --> " << __LINE__ << std::endl;
      ofs << "# File : " << fileCreatedIn << std::endl;
      ofs << "# Function : " << functionCreatedIn << std::endl;
      ofs << "# Line Number : " << lineNumber << std::endl;
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __FUNCTION__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR");
      exit(EXIT_FAILURE);
    }
  
  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
void makeInfo(std::vector<std::string>& info, const std::string& fileName,
	      const char* const file, const char* const func, const unsigned& line)
{
  info.resize(4);
  info[0] = fileName;
  info[1] = file;
  info[2] = func;
  info[3] = num2str(line);
  return;
}
//------------------------------






/*
//--------------------------------------------------------------------------------//
void sortRows(double* const array, const unsigned & numRows, const unsigned & numCols, 
	      const unsigned & keyCol)
{
  if ( keyCol > numCols-1 )
    {
      cerr << "File: " << __FILE__ << endl;
      cerr << "Function: " << __PRETTY_FUNCTION__  << endl;
      cerr << "ERROR: " 
	   << "key column can not be greater than the number of columns"
	   << endl;
      exit(EXIT_FAILURE);
    }

  double* tempArray = makeArray(tempArray, numCols);


  double num1 = 0.0;
  double num2 = 0.0;


  unsigned minValueRow = 0; 
  // selection sort
  for( unsigned n = 0; n < (numRows-1); ++n)
    {
      // assume that the smallest
      // element in the key column is at row
      minValueRow = n;

      // search over remaining elements to find
      // smallest element
      for(unsigned m = n+1; m < numRows; ++m)
	{
	  // check if a smaller value exists
	  num1 = array[m*numCols+keyCol];
	  num2 = array[minValueRow*numCols+keyCol];
	  if ( gsl_fcmp( num1, num2, DBL_EPSILON ) == -1 )
	    {
	      // save index value
	      // for later swap
	      minValueRow = m;
	    }
	}

      // swap current with that at minValueRow

      //tempInd = Ind[minValueIndex];
      for( unsigned col = 0; col < numCols; ++col)
	{
	  tempArray[col] = array[minValueRow*numCols+col];
	}


      //Ind[minValueIndex] = Ind[n];
      for( unsigned col = 0; col < numCols; ++col)
	{
	  array[minValueRow*numCols+col] = array[n*numCols+col];
	}

      //Ind[n] = tempInd;
      for( unsigned col = 0; col < numCols; ++col)
	{
	  array[n*numCols+col] = tempArray[col];
	}
    }


  deleteArray(tempArray);
  return;
}
//--------------------------------------------------------------------------------//
*/

//--------------------------------------------------------------------------------//
void sortRows(double* const array, size_t const numRows, size_t const numCols, 
	      size_t const keyCol)
{
  assert(keyCol < numCols);

  double* tempArray = makeArray<double>(numRows*numCols);

  size_t* p = makeArray<size_t>(numRows);
  const unsigned stride = numCols;
  gsl_sort_index(p, array+keyCol, stride, numRows);

  for(unsigned row = 0; row < numRows; ++row)
    {
      unsigned index = p[row];
      for(unsigned col = 0; col < numCols; ++col)
	{
	  tempArray[row*numCols+col] = array[index*numCols+col];
	}
    }

  for(unsigned n = 0; n < numRows*numCols; ++n)
    {
      array[n] = tempArray[n];
    }

  
  delete [] tempArray;
  delete [] p;
  return;
}
//------------------------------



//--------------------------------------------------------------------------------//
void dec2bin(const double & xmax, const double & xmin, const double & x, double & approx, 
	     unsigned* const bitArray, const unsigned & numBits)
{
  // x is the number that needs to be approximated
  // approx is the approximation of x

  assert(xmin < xmax);

  const double twoPower = pow(2.0, double(numBits));

  // how many levels or codes
  const double numLevels = (twoPower-1);

  // resolution
  // difference between 2 consecutive representations
  const double delta  = (xmax - xmin)/numLevels;

  const double deltaByTwo = delta/2.0;

  if((x < (xmin-deltaByTwo)) || (x > (xmax+deltaByTwo)))
    {
      std::cerr << __FILE__ << " --> " << __FUNCTION__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR :: ");
      std::cerr << "x not in between (xmax  xmin)" << std::endl; 
      std::cerr << "x = " << x << std::endl;
      std::cerr << "xmax = " << xmax << std::endl;
      std::cerr << "xmin = " << xmin << std::endl;
      exit(EXIT_FAILURE);
    }

  // which level does x belong to ?
  double howmanydeltas = (x-xmin)/delta;

  unsigned howMany = floor(howmanydeltas+0.5);
  // floor(x+0.5) == round(x) 
  // http://faculty.salisbury.edu/~dxdefino/RoundOff.htm
  // there are many ways of rounding
  // http://en.wikipedia.org/wiki/Rounding
  // http://stackoverflow.com/questions/485525/round-for-float-in-c

  approx = xmin + howMany*(xmax-xmin)/numLevels;

  for(unsigned i = 0; i < numBits; ++i)
    {
      // AND with 1 and extract the lsb of the product
      bitArray[i] =  ( (howMany >> i) & 0x1 );
    }

  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
void dec2bin(bool & result, const double & xmax, const double & xmin, const double & x, double & approx, 
	     unsigned* const bitArray, const unsigned & numBits)
{
  // x is the number that needs to be approximated
  // approx is the approximation of x

  assert(xmin < xmax);

  const double twoPower = pow(2.0, double(numBits));

  // how many levels or codes
  const double numLevels = (twoPower-1);

  // resolution
  // difference between 2 consecutive representations
  const double delta  = (xmax - xmin)/numLevels;

  const double deltaByTwo = delta/2.0;

  if((x < (xmin-deltaByTwo)) || (x > (xmax+deltaByTwo)))
    {
      result = false;
    }
  else
    {
      result = true;
      
      // which level does x belong to ?
      double howmanydeltas = (x-xmin)/delta;
      
      unsigned howMany = floor(howmanydeltas+0.5);
      // floor(x+0.5) == round(x) 
      // http://faculty.salisbury.edu/~dxdefino/RoundOff.htm
      // there are many ways of rounding
      // http://en.wikipedia.org/wiki/Rounding
      // http://stackoverflow.com/questions/485525/round-for-float-in-c
      
      approx = xmin + howMany*(xmax-xmin)/numLevels;
      
      for(unsigned i = 0; i < numBits; ++i)
	{
	  // AND with 1 and extract the lsb of the product
	  bitArray[i] =  ( (howMany >> i) & 0x1 );
	}
    }
  return; 
}
//------------------------------



//--------------------------------------------------------------------------------//
void bin2dec(const unsigned* const bitArray, const unsigned & len,
	     unsigned & decimalNumber)
{
  decimalNumber = 0;

  int n = len-1;

  while(n >= 0)
    {
      // 2*x == x<<1 
      // left shift by one = multiplication by 2
      decimalNumber = decimalNumber<<1;
      decimalNumber = decimalNumber + bitArray[n];
      n = n-1;
    }

  return;
}
//--------------------------------------------------------------------------------//


//--------------------------------------------------------------------------------//
unsigned bin2dec(const unsigned* const bitArray, const unsigned & len)
{
  unsigned decimalNumber = 0;

  int n = len-1;

  while(n >= 0)
    {
      // 2*x == x<<1 
      // left shift by one = multiplication by 2
      decimalNumber = decimalNumber<<1;
      decimalNumber = decimalNumber + bitArray[n];
      n = n-1;
    }

  return decimalNumber;
}
//--------------------------------------------------------------------------------//


//--------------------------------------------------------------------------------//
bool isEqual(const double * const array1, const double * const array2, 
	     size_t const size){
  bool result = true;
  for(size_t element = 0; element < size; ++element){
    if(gsl_fcmp( array1[element], array2[element], DBL_EPSILON ) != 0){
      result = false;
      break;
    }
  }
  return result;
}
//------------------------------


//------------------------------------------------------------------------------------------
// check if two arrays of real numbers are approximtely equal
bool isEqual(const std::vector<double>& array1, const std::vector<double>& array2)
{
  assert(array1.size()==array2.size()); //?????
  
  bool result = true;
  for(unsigned element = 0; element < array1.size(); ++element)
    {
      if(gsl_fcmp(array1[element], array2[element], DBL_EPSILON) != 0) 
	{
	  result = false;
	  break;
	}
    }

  return result;
}
//------------------------------

//------------------------------------------------------------------------------------------
//Sat Jul 31 08:09:50 2010
bool equal(const double & x, const double & y)
{
  if(gsl_fcmp(x, y, DBL_EPSILON) == 0) {return true;}else{return false;}
}
//------------------------------


//------------------------------------------------------------------------------------------
//Sat Aug  7 01:09:48 2010
bool lessThan(const double & x, const double & y)
{
  if(gsl_fcmp(x, y, DBL_EPSILON) == -1) 
    {// x < y
      return true;
    }
  else
    {// y < x
      return false;
    }
}
//------------------------------

//------------------------------------------------------------------------------------------
//Sat Aug  7 06:53:46 2010
bool greaterThan(const double & x, const double & y)
{
  if(gsl_fcmp(x, y, DBL_EPSILON) == 1) 
    {// x > y
      return true;
    }
  else
    {// y > x
      return false;
    }
}
//------------------------------


//------------------------------------------------------------------------------------------
bool fileExists(const std::string & fileName)
{
  // check if a file exists 
  // http://www.techbytes.ca/techbyte103.html
  
  struct stat fileStats;
  int intStat = stat(fileName.c_str(), &fileStats);

  //bool returnValue = false;

  if(intStat == 0)
    {
      // file exists
      //returnValue = true; 
      std::cerr << alert("ERROR:: ") << "file " << fileName << " already exists ... " << std::endl;
    }
  else
    {
      // file does not exist
      //returnValue = false; 
    }

  //return returnValue;
  return false; ////// WRONG  WRONG WRONG WRONG WRONG WRONG WRONG 
}
//------------------------------

//------------------------------------------------------------------------------------------
void generateIndexLists(const gsl_rng* grng, unsigned* & list1, unsigned* & list2, 
			const unsigned & size)
{
  // with list length 2 can not avoid ab ba situation 
  assert(size != 2);

  list1 = new unsigned [size];
  list2 = new unsigned [size];

  for(unsigned i = 0; i < size; ++i){list1[i] = list2[i] = i;}

  gsl_ran_shuffle(grng, list1, size, sizeof(unsigned));
  gsl_ran_shuffle(grng, list2, size, sizeof(unsigned));

  // make sure numbers occupying identical places in two lists
  // are not identical
  bool condition1 = true;
  
  // also check if we have a (5,2) and (2,5) situation i.e. both 2 and 5
  // participate in exactly two tournaments but each time they fight
  // against each other ... for the sake of diversity and variability we
  // do not want this
  bool condition2 = true; 

  do
    {
      condition1 = true;
      for(unsigned i = 0; i < size; ++i)
	{
	  if(list1[i] == list2[i])
	    {
	      condition1 = false; break;
	    }
	}

      if(condition1)
	{
	  condition2 = true;
	  std::vector<bool> ignore(size, false);
	  for(unsigned i = 0; i < size; ++i)
	    {
	      if(ignore[i]){continue;}
	      
	      for(unsigned j = i+1; j < size; ++j)
		{
		  if(ignore[j]){continue;}
		  
		  if(list1[i] == list2[j] && list2[i] == list1[j])
		    {
		      condition2 = false;
		      ignore[j] = true;
		    }
		}
	    }
	}

      if(!condition1 || !condition2)
	{
	  gsl_ran_shuffle(grng, list2, size, sizeof(unsigned));
	}
      
    }while(!(condition1 && condition2));

  return;
}
//----------------------------------------------------------------------


//------------------------------------------------------------------------------------------
std::string alert(const std::string & str)
{
  return "\033[0;31m" + str + "\033[0m";
}
//------------------------------
// not entirely happy with the implementation here ... first a temporary
//object is created and then i m passing it back by copy ... :( not good
//----------------------------------------------------------------------


//------------------------------------------------------------------------------------------
std::string ingreen(const std::string & str)
{
  return greenColor + str + noColor;
}
//------------------------------

//------------------------------------------------------------------------------------------
std::string inred(const std::string & str)
{
  return redColor + str + noColor;
}
//------------------------------


//----------------------------------------------------------------------
void warning()
{
  //std::cout << "\a";
  std::cout << "\033[0;31mwarning::\033[0m" << std::endl;
  return;
}
//----------------------------------------------------------------------

//----------------------------------------------------------------------
// default argument is set in the header file :)
void newline(const size_t numLines){
  for(size_t n = 0; n < numLines; ++n)
    {
      std::cout << "\n";
    }
  return;
}
//----------------------------------------------------------------------


//------------------------------------------------------------------------------------------
// void paretoSort(const std::vector<std::vector<double> >& data,
// 		std::vector<unsigned>& best)
// {
//   // takes the candidates for gapFillers and finds out the best
//   // non-dominated set amongst them
  
//   const unsigned numVec = data.size();
  
//   if(numVec != best.size())
//     {
//       std::cerr << "File: " << __FILE__ << std::endl;
//       std::cerr << "Function: " << __FUNCTION__  << std::endl;
//       std::cerr << "Line: " << __LINE__  << std::endl;
//       std::cerr << alert("ERROR:: ") 
// 		<< "size mismatch " << std::endl;   
//       exit(EXIT_FAILURE);
//     }
  
//   for(unsigned n = 0; n < numVec; ++n)
//     {
//       best[n] = n;
//     }
  
//   std::vector<unsigned>::iterator iter;
  
//   for(iter = best.begin(); iter != best.end(); ++iter)
//     {
//       const unsigned one = *iter;

//       for(unsigned two = (one+1); two < numVec; ++two)
// 	{
// 	  const int result = dominanceCompare(data[one], data[two]);
	  
// 	  if(result==-1)
// 	    {
// 	      // two dominates one
// 	      // remove one from best
// 	      best.erase(iter);
// 	      //erase invalidates all iterator and references to elements after iter.
// 	      iter--;
// 	      break;
// 	    }
// 	  else if(!result)
// 	    {
// 	      // both non-dominated
// 	    }
// 	  else if(result==1)
// 	    {
// 	      // one dominates two
// 	      // remove two from best
// 	      best.erase(std::remove(best.begin(), best.end(), two), best.end());
// 	    }
// 	  else
// 	    {
// 	      std::cerr << "this should not happen" << std::endl;
// 	      break;
// 	    }
// 	}
//     }

//   return;
// }
//------------------------------
