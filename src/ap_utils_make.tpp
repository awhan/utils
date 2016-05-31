#ifndef __ap_utils_make_tpp__
#define __ap_utils_make_tpp__

#include <iostream>

//--------------------------------------------------------------------------------//
template <typename T>
T *** makeArray(size_t const numTables, size_t const numRows, size_t const numCols)
{
  T *** array = new double** [numTables];
  if ( array == NULL )
    {
      std::cerr << "File: " << __FILE__ << std::endl;
      std::cerr << "Function: " << __FUNCTION__ << std::endl;
      std::cerr << "Error: can't allocate memory for "
		<< "array" << std::endl;
      exit(EXIT_FAILURE);
    }
  else
    {
      for(unsigned table = 0; table < numTables; ++table)
	{
	  array[table] = new double*[numRows];
	  if ( array[table] == NULL )
	    {
	      std::cerr << "File: " << __FILE__ << std::endl;
	      std::cerr << "Function: " << __FUNCTION__ << std::endl;
	      std::cerr << "Error: can't allocate memory for " 
			<< "array[" << table << "]" << std::endl;
	      exit(EXIT_FAILURE);
	    }
	  else
	    {
	      for(unsigned row = 0; row < numRows; ++row)
		{
		  array[table][row] = new double [numCols];
		  if ( array[table][row] == NULL )
		    {
		      std::cerr << "File: " << __FILE__ << std::endl;
		      std::cerr << "Function: " << __FUNCTION__ << std::endl;
		      std::cerr << "Error: can't allocate memory for "
			   << "array[" << table << "][" 
				<< row << "]" << std::endl;
		      exit(EXIT_FAILURE);
		    }
		  else 
		    {
		      // should i initialize to zero ???
		      for(unsigned col = 0; col < numCols; ++col)
			{
			  array[table][row][col] = 0.0;
			}
		    }
		}
	    }
	}
    }


  return array;
}
//--------------------------------------------------------------------------------//

template <typename T>
void deleteArray(T *** array, size_t const numTables, size_t const numRows) {
  for(size_t table = 0; table < numTables; ++table) {
      for(size_t row = 0; row < numRows; ++row) {
        delete [] array[table][row];
      }
      delete [] array[table];
  }
  delete [] array;
}
//--------------------------------------------------------------------------------//


//--------------------------------------------------------------------------------
template <class T>
T * makeArray(size_t const  numElements)
{
  T * array = new T [numElements];
  if(!array)
    {
      std::cerr << "ERROR:: " << "can't allocate memory ..." << std::endl;
      exit(EXIT_FAILURE);
    }
  
  return array;
}
//------------------------------

//------------------------------------------------------------------------------------------
template <class t>
t** makeArray(const unsigned & numRows, const unsigned & numCols)
{
  t** array = new t* [numRows];

  if(!array)
    {
      //std::cerr << alert("ERROR:: ") << "can't allocate memory ..." << std::endl;
      std::cerr << "ERROR:: " << "can't allocate memory ..." << std::endl;
      exit(EXIT_FAILURE);
    }
  
  for(unsigned rowIndex = 0; rowIndex < numRows; ++rowIndex)
    {
      array[rowIndex] = new t [numCols];
      
      if(!(array+rowIndex))
	{
	  std::cerr << "ERROR:: " << "can't allocate memory ..." << std::endl;
	  exit(EXIT_FAILURE);
	}
    }

  return array;
}
//--------------------------------------------------------------------------------/

//--------------------------------------------------------------------------------//
template <typename T>
void deleteArray(T ** array, size_t const numRows) {
  for(size_t row = 0; row < numRows; ++row) {
      delete [] array[row];
    }
  delete [] array;
}
//--------------------------------------------------------------------------------//

#endif
