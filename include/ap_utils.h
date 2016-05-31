#ifndef __ap_utils_h__
#define __ap_utils_h__

//

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h> // for axblu();
#include <gsl/gsl_cblas.h>

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iosfwd>

#include <cassert>


//#include "ap_utils_gsl.h"


///////////////////////////////////////////////////////////////////////////////

#define SHOW(x) std::cout << #x " :: " << x << std::endl


// #define SHOW(x) std::cout << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n" << #x " :: " << x << std::endl

//#define apSHOW(x, ...) std::cout << #x << std::endl; display(x
//#,##__VA_ARGS__)

//#define ap_show(x) std::cout << "\n" __FILE__ << " --> " << __func__ << "()
//#--> " << __LINE__ << "\n"; display(x);

#define AP_SHOW(x) std::cout << #x << std::endl; display(x)
#define showv(x) std::cout << #x << "\t"; display(x)


const std::string noColor	= "\033[0m";
std::string const bold		= "\033[1m";

const std::string redColor	= "\033[0;31m";
const std::string greenColor	= "\033[0;32m";
const std::string yellowColor	= "\033[0;33m";
const std::string blueColor	= "\033[0;34m";
const std::string magentaColor	= "\033[0;35m";
const std::string cyanColor	= "\033[0;36m";

std::string const boldred       = "\033[1;31m";
std::string const boldgreen     = "\033[1;32m";
std::string const boldyellow    = "\033[1;33m";
std::string const boldblue      = "\033[1;34m";
std::string const boldmagenta   = "\033[1;35m";
std::string const boldcyan      = "\033[1;36m";
std::string const boldwhite     = "\033[1;37m";


#if 0
///////////////////////////////////////////////////////////////////////////////
template <typename T> std::vector<T> read(std::string const filename, size_t
                                          const nrows, size_t const ncols) {
  std::cout << "\n" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
  std::ifstream ifs(filename.c_str(), std::ios::in);
  std::cout << nrows << ncols << std::endl;
  std::string line;
  while(std::getline(ifs, line)) {
    std::cout << line << std::endl;
    std::istringstream iss(line);
    std::string w;
    while(iss >> w) { std::cout << w << std::endl; }
  }
  std::vector<T> v;
  return v;
}
#endif


//------------------------------------------------------------------------------------------
template <typename t> void
transpose(t* out, t* const in, size_t const nrows, size_t const ncols)
{
  for(size_t i = 0; i < ncols; ++i)
    for(size_t j = 0; j < nrows; ++j)
      out[i*nrows+j] = in[j*ncols+i];
}
//------------------------------

//------------------------------------------------------------------------------------------
template <typename t> void
transpose(t* matrix, size_t const nrows, size_t const ncols)
{
  size_t const numel = nrows*ncols;
  t* temp = new t[numel];
  
  for(size_t i = 0; i < ncols; ++i)
    for(size_t j = 0; j < nrows; ++j)
      temp[i*nrows+j] = matrix[j*ncols+i];
  
  std::copy(temp, temp+numel, matrix);
  
  delete [] temp;
}
//------------------------------


//size_t rank(gsl_matrix * const m, double const tol = 0.0);
size_t rank(gsl_matrix * const m);

std::string ingreen(const std::string & str);
std::string inred(const std::string & str);


// analogous to gnu octave center()
void center(std::vector<std::vector<double> >& matrix, const int flag);
void center(gsl_matrix* matrix, gsl_vector* const means);
  
// analogous to gnu octave mean()
void mean(std::vector<double>& means, const std::vector<std::vector<double> >& matrix,
	  const int flag);
void mean(gsl_vector* means, const gsl_matrix* matrix);

// compute the euclidean distance matrix 
void distanceMatrix(std::vector<std::vector<double> >& dm,
		    const std::vector<std::vector<double> >& data,
		    const std::vector<std::vector<double> >& centers);

void distanceMatrix(double* out, double* const data, size_t const numData, double* const centers,
		    size_t const numCenters, size_t const ndims);

bool nonDominated(const std::vector<double>& v, const std::vector<std::vector<double> >& a);


// compute the pseudo inverse of a matrix, square, fat or thin
void pinv(gsl_matrix* apinv, const gsl_matrix* a);

// solve for X in AX=B by X = pinv(A)*B method
void axbpinv(const std::vector<std::vector<double> >& a,
	     std::vector<std::vector<double> >& x,
	     const std::vector<std::vector<double> >& b,
	     const unsigned& m, const unsigned& n, const unsigned& p);


// solve for X in AX=B by LU decomposition method. 
void axblu(const std::vector<std::vector<double> >& a, std::vector<std::vector<double> >& x,
	   const std::vector<std::vector<double> >& b, size_t const m, size_t const p);

// convert STLGSL to GSL 
void stl2gsl(const std::vector<std::vector<double> >& stl, gsl_matrix* gsl,
	     const unsigned& numRows, const unsigned& numCols);

// convert GSL to STL 
void gsl2stl(const gsl_matrix* gsl, std::vector<std::vector<double> >& stl,
	     const unsigned& numRows, const unsigned& numCols);


template <class t>
std::vector<t> subtract(const std::vector<t>& a, const std::vector<t>& b);


void makeFile(std::ofstream& ofs, const std::vector<std::string>& info);

void makeInfo(std::vector<std::string>& info, const std::string& fileName,
	      const char* const file, const char* const func, const unsigned& line);

// template <typename t>
// void display(const t& data);

// single value 
template <typename t>
void write(t const value, const std::string& fileName);

// 4-dim
template <typename t>
void write(const std::vector<std::vector<std::vector<std::vector<t> > > >& array, 
	   const std::vector<std::string>& info, const unsigned& togetherLines = 1); 

// 3-dim
template <typename t>
void write(const std::vector<std::vector<std::vector<t> > >& array, 
	   const std::vector<std::string>& info);

// 3-dim to file stream
template <typename t>
void write(const std::vector<std::vector<std::vector<t> > >& array, std::ofstream& ofs,
	   const std::vector<std::string>& info);

template <class t>
void write(const std::vector<t>& array, std::ofstream& ofs); 

template <class t>
void write(const std::vector<std::vector<t> >& array, const std::string& fileName);

template <class t>
void write(const std::vector<std::vector<t> >& array, std::ofstream& ofs,
	   const std::vector<std::string>& info);

template <class t>
void write(const t* const data, const unsigned& numRows, const unsigned& numCols,
	   const std::vector<std::string>& info);

template <class t>
void write(const std::vector<std::vector<t> >& array, const std::vector<std::string>& info);

// 2-dim
template <typename t>
void display(const std::vector<std::vector<t> >& array);

// 3-dim
template <typename t>
void display(const std::vector<std::vector<std::vector<t> > >& array);

// 1-dim and make a value stand out in the list
template <typename t>
void display(const t& array, const typename t::value_type& value);

template <typename t>
void display(const t* const array, size_t const numElements);

template <typename t>
void display(const t* const array, const unsigned numRows, const unsigned numCols);

template <typename t>
void display(const std::vector<t>& array, std::size_t const from, std::size_t const to);

//++++++++++ 1 Dimensional Arrays++++++++++//


// quick and dirty testing
template <typename t>
void write(const std::vector<t>& array, const std::string& fileName);

template <typename T> T * makeArray(size_t const numElements);

template <class t>
t** makeArray(const unsigned& numRows, const unsigned& numCols);

// 2-dim array 
// this is used for quick and dirty testing
//Sun Sep  5 22:36:37 2010
template <class t>
void write(const t* const array, const size_t size, const size_t stride, const size_t numRows,
	   const size_t numCols, const std::string& fileName);

void read(const std::string& fileName, double* data, const unsigned numRows, const unsigned numCols);


// read an array from a file
void read(const std::string& fileName, std::vector<std::vector<double> >& data, 
	  const unsigned& numRows, const unsigned& numCols);


//------------------------------------------------------------------------------------------
template <typename T>
std::string num2str(T const num)
{
  std::stringstream ss;
  ss << num; 
  return ss.str();
}
//------------------------------


////////////////////////////////////////////////////////////////////////////////
template <typename T>
T str2num(std::string const str) {
  std::istringstream iss(str);
  T num{};
  iss >> num;
  return num;
}

// sort the rows of a matrix according to any column
void sortRows(double* const array, const unsigned& numRows, const unsigned& numCols, 
	      const unsigned& keyCol);

template <typename t>
void sortRows(std::vector<std::vector<t> >& array, const unsigned& keyColumn);

//Sat Jul 24 00:48:06 2010
// generate a sorted index list of matrix rows by a specified key column
template <typename t>
void sortIndex(std::vector<unsigned>& srtIdx, const std::vector<std::vector<t> >& array,
	       const unsigned& keyColumn);


void dec2bin(bool& result, const double& xmax, const double& xmin, const double& x, double& approx, 
	     unsigned* const bitArray, const unsigned& numBits);

void dec2bin(const double& xmax, const double& xmin, 
	     const double& x, double& approx, unsigned* const bitArray,
	     const unsigned& len);

void bin2dec(const unsigned* const bitArray, const unsigned& len,
	     unsigned& decimalNumber);

unsigned bin2dec(const unsigned* const bitArray, const unsigned& len);

// check if two arrays of real numbers are approximtely equal
bool isEqual(const double* const array1, const double* const array2, 
	     const unsigned& size);

// check if two arrays of real numbers are approximtely equal
bool isEqual(const std::vector<double>& array1, const std::vector<double>& array2);


//Sat Jul 31 08:09:00 2010
bool equal(const double&, const double&);

//Sat Aug  7 01:08:24 2010
bool lessThan(const double&, const double&);

//Sat Aug  7 06:53:00 2010
bool greaterThan(const double&, const double&);


// check if a file exists
bool fileExists(const std::string& fileName);

// check dominance relation between two vectors
// A dominates B :: return 1
// B dominates A :: return -1
// A = B or A and B are non-dominant :: return 0
int dominanceCompare(const double* const A, const double* const B,
		     const unsigned& size);

// A dominates B :: return 1
// B dominates A :: return -1
// A = B or A and B are non-dominant :: return 0
int dominanceCompare(const std::vector<double>& A, const std::vector<double>& B);

void generateIndexLists(const gsl_rng*, unsigned*& list1, 
			unsigned*& list2, const unsigned& size);

std::string alert(const std::string& str);

// error is a standard name in glibc ... consider changing it !!!
std::string error();
void warning();


// default argument value 1
void newline(size_t const numLines = 1);

void paretoSort(const std::vector<std::vector<double> >& data,
		std::vector<unsigned>& best);

// to count the number of vectors in a 3-dim array, if you have a
// vector<vector<vector<double> > > d3 then you can count the number
// of vectors it has. this is a function object aka functor to be used
// in conjunction with algorithms like for_each
template <class T>
class vectorCounter
{
public:
  // constructor
  vectorCounter(unsigned& numElements):numElements_(numElements){}
  
  void operator()(const T& iter)
  {
    numElements_ = numElements_ + iter.size();
  }
  
private:
  unsigned& numElements_;
};


#include "ap_utils_math.tpp"
#include "ap_utils_display.tpp"
#include "ap_utils_sort.tpp"
#include "ap_utils_make.tpp"
#include "ap_utils_write.tpp"
#include "ap_utils_stats.tpp"
#include "ap_utils_read.tpp"


#endif
