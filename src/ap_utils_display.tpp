#ifndef ap_utils_display_tpp__
#define ap_utils_display_tpp__

#include <iostream>

//#define AP_VERBOSE_DISP

//------------------------------------------------------------------------------------------
template <typename T>
void display(const T & arg) {
  #ifdef AP_VERBOSE_DISP
  std::cout << "# " << arg.size() << " elements " << " :: " << __FILE__ << " --> "
	    << __func__ << " --> " << __LINE__ << "\n";
  #endif
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  for(auto & i : arg) {
    std::cout << i << "\t";
  }
  std::cout << "\n";
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
}
//------------------------------

//------------------------------------------------------------------------------------------
template <typename InputIterator>
void display(InputIterator begin, InputIterator end) {
  //std::cout << "# " << distance(begin, end) << " elements\n";
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  while(begin != end) {
    std::cout << *begin++ << "\t";
  }
  std::cout << std::endl;
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
}
//------------------------------  

//------------------------------------------------------------------------------------------
template <typename t>
void display(const std::vector<std::vector<t> >& array) {
  #ifdef AP_VERBOSE_DISP
  std::cout << "\n" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
  #endif

  typedef typename std::vector<std::vector<t> >::size_type sizeType1;
  //typedef typename decltype(array)::size_type sizeType1;
  const sizeType1 numRows = array.size();
  typedef typename std::vector<t>::size_type sizeType2;
  
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  std::cout << "# " << numRows << " x " << array[0].size() << std::endl;
  for(sizeType1 row = 0; row < numRows; ++row)
    {
      // row number
      std::cout << row << ")\t";

      for(sizeType2 col = 0; col < array[row].size(); ++col)
	{
	  if(col == (array[row].size()-1))
	    {
	      std::cout << array[row][col];
	    }
	  else 
	    {
	      std::cout << array[row][col] << "\t";
	    }
	}
      std::cout << std::endl;
    }
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint); 
}
//------------------------------

//------------------------------------------------------------------------------------------
template <typename t>
void display(const std::vector<std::vector<std::vector<t> > >& array)
{
  std::cout << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  for(unsigned table = 0; table < array.size(); ++table)
    {
      std::cout << "# table " << table << std::endl;
      
      for(unsigned row = 0; row < array[table].size(); ++row)
	{
	  std::cout << row << ")\t";
	  
	  for(unsigned col = 0; col < array[table][row].size(); ++col)
	    {
	      if(col == (array[table][row].size()-1))
		{
		  std::cout << array [table][row][col];
		}
	      else 
		{
		  std::cout << array [table][row][col] << "\t";
		}
	    }
	  std::cout << std::endl;
	}
      if(table == (array.size()-1))
	{}
      else
	{
	  std::cout << std::endl;
	}
    }
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  return;
}
//--------------------------------------------------


template <typename T> 
void display(T *** array, size_t const npages, size_t const nrows, size_t const ncols) {
  std::cout << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  for(size_t p = 0; p < npages; ++p) {
    for(size_t r = 0; r < nrows; ++r) {
      for(size_t c = 0; c < ncols; ++c) {
        std::cout << array[p][r][c] << "\t"; // FIXME :: no trailing TABs
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  return;
}


//------------------------------------------------------------------------------------------
template <typename t> void
display(const t * const array, size_t const numElements){
  // display array
  assert(array);
  std::cout << "# " << numElements << " elements " << " :: " << __FILE__ << " --> "
	    << __func__ << " --> " << __LINE__ << "\n";
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  for(size_t index = 0; index < numElements; ++index){
    if(index == (numElements-1)) {
      std::cout << array[index] << std::endl;
    } else {
      std::cout << array[index] << "\t";
    }
  }
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint); 
}
//------------------------------

// //------------------------------------------------------------------------------------------
// template <typename t> 
// void display(const t * const array, size_t const numRows, size_t const numCols){
//   assert(array);
//   std::cout << "# " << numRows << " x " << numCols << " :: " << __FILE__ << " --> "
// 	    << __func__ << "() --> " << __LINE__ << "\n";
//   std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
//   for(size_t row = 0; row < numRows; ++row){
//     std::cout << row << ")\t";
//     for(size_t col = 0; col < numCols; ++col){
//       if(col == (numCols-1)){
//         std::cout << array[row*numCols+col] << std::endl;
//       }else{
//         std::cout << array[row*numCols+col] << "\t";
//       }
//     }
//   }
//   std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint); 
// }
// //------------------------------



//------------------------------------------------------------------------------------------
template <typename in_itr>
void display(in_itr ii, size_t const nrows, size_t const ncols){
  std::cout << "# " << nrows << " x " << ncols << " :: " << __FILE__ << " --> "
	    << __func__ << "() --> " << __LINE__ << "\n";
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  for(size_t row = 0; row < nrows; ++row){
    //std::cout << row << ")\t";
    for(size_t col = 0; col < ncols; ++col){
      if(col == (ncols-1)){
        std::cout << *(ii++) << std::endl;
      }else{
        std::cout << *(ii++) << "\t";
      }
    }
  }
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint); 
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t>
void display(const t& array, const typename t::value_type& value)
{
  // make a value stand out in the list
  
  const unsigned size = array.size();
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  std::cout << "# " << size << " elements " << std::endl;
  const auto last = std::prev(array.cend());
  auto itr = array.cbegin();
  for(; itr != last; ++itr)
    {
      if(value == *itr)
	std::cout << redColor << *itr << noColor << ", ";
      else
	std::cout << *itr << ",  ";
    }
  if(value == *itr)
    std::cout << redColor << *itr << noColor << std::endl;
  else
    std::cout << *itr << std::endl;
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint); 

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t>
void display(const std::vector<t>& array, std::size_t const from, std::size_t const to)
{
  assert(from < to);
  assert(to < array.size());
  
  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint); 
  for(unsigned index = from; index < to; ++index)
    {
      if(index == (array.size()-1))
	{
	  // to endl or not to endl ??? 
	  std::cout << array[index] << std::endl;
	}
      else 
	std::cout << array[index] << "\t";
    }
  std::cout.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint); 

  return;
}
//------------------------------

#endif
