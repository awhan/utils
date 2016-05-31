#ifndef _SR_H_
#define _SR_H_

#include <vector>
#include <algorithm>

#include <gsl/gsl_sort.h>
// for gsl_sort_index()

template <typename t>
struct comp
{
  typename std::vector<t>::size_type idx;
  comp(const typename std::vector<t>::size_type & col) : idx(col) { }
  
  bool operator()(std::vector<t> const& a, std::vector<t> const& b)
  {
    return a.at(idx) < b.at(idx);
  }
};


template <typename t>
void sortRows(std::vector <std::vector <t> > & array, const unsigned int & keyColumn)
{
  std::sort(array.begin(), array.end(), comp<t>(keyColumn));
  return;
}

// TODO
// template <typename t>
// struct compIndex
// {
//   std::vector <t> arr;
//   compIndex(const std::vector <t> & array): arr(array){}

//   bool operator()
// }

// template <typename t>
// void sortIndex(const std::vector <t> & array, std::vector <unsigned int> & indexList)
// {
//   std::sort(array.begin(), array.end(), compIndex<t>(keyColumn));
//   return;
// }


//------------------------------------------------------------------------------------------
//Sat Jul 24 00:48:19 2010
template <typename t>
void sortIndex(std::vector <unsigned> & srtIdx, const std::vector <std::vector <t> > & array,
	       const unsigned & keyColumn)
{
  const unsigned numRows = array.size();

  t * data = makeArray<t>(numRows);
  
  for(unsigned row = 0; row < numRows; ++row)
    { 
      data[row] = array[row][keyColumn];
    }

   unsigned * srtIdxTemp = makeArray<unsigned>(numRows);
   const unsigned stride = 1;
   gsl_sort_index(srtIdxTemp, data, stride, numRows);
   delete [] data;

   std::copy(srtIdxTemp, srtIdxTemp+numRows, srtIdx.begin());

   delete [] srtIdxTemp;
   return;
}
//------------------------------

#endif 
