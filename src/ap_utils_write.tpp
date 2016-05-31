#ifndef ap_utils_write_tpp__
#define ap_utils_write_tpp__

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <ios>


//------------------------------------------------------------------------------------------
template <typename t>
void write(const std::vector<std::vector<std::vector<std::vector<t> > > >& array,
           const std::vector<std::string>& info, const unsigned& togetherLines)
{
  // 4-dimensional array[dim0][dim1][dim2][dim3]

  const std::string fileName = info[0];
  const std::string fileCreatedIn = info[1];
  const std::string functionCreatedIn = info[2];
  const std::string lineNumber = info[3];

  assert(!fileExists(fileName));

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  //const unsigned togetherLines = 2;

  if(outFile.good())
    {
      outFile << "# " << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      outFile << "# File : " << fileCreatedIn << "\n";
      outFile << "# Function : " << functionCreatedIn << "\n";
      outFile << "# Line : " << lineNumber << "\n";
      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      outFile.precision(std::numeric_limits<t>::digits10);

      for(unsigned n = 0, idx0 = 0; idx0 < array.size(); ++idx0)
        {
          outFile << "\n# dim0 = " << idx0 << "\n";
          for(unsigned idx1 = 0; idx1 < array[idx0].size(); ++idx1)
            {//table
              outFile << "## table = " << idx1 << " of dim0 = " << idx0 << "\n";
              for(unsigned idx2 = 0; idx2 < array[idx0][idx1].size(); ++idx2, ++n)
                {//row
                  outFile << n << "\t";
                  for(unsigned idx3 = 0; idx3 < array[idx0][idx1][idx2].size(); ++idx3)
                    {//col
                      if(idx3 == (array[idx0][idx1][idx2].size()-1))
                        {
                          outFile << array[idx0][idx1][idx2][idx3] << "\n";
                        }
                      else
                        {
                          outFile << array[idx0][idx1][idx2][idx3] << "\t";
                        }
                    }
                  // double space
                  if((idx2+1)%togetherLines==0)
                    {
                      if(idx2 !=(array.size()-1)){outFile << "\n\n";}
                    }
                }
              if(idx1 != array[idx0].size()-1)outFile << "\n";
            }
        }
      outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      std::cerr << alert("ERROR: ")
                << "could not open file " << fileName << std::endl;
      exit(EXIT_FAILURE);
    }

  outFile.close();

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t>
void write(const std::vector<std::vector<t> >& array, const std::vector<std::string>& info,
           const unsigned& togetherLines)
{
  // the following will write the array but introduce double lines after
  // writing togetherLines number or lines if togetherLines=4 then there
  // will be a double space after every 4 lines, this is needed to use
  // the `using index' feature of gnuplot

  const std::string fileName = info[0];
  const std::string fileCreatedIn = info[1];
  const std::string functionCreatedIn = info[2];
  const std::string lineNumber = info[3];

  assert(!fileExists(fileName));

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good())
    {
      outFile << "#" << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      outFile << "# File : " << fileCreatedIn << "\n";
      outFile << "# Function : " << functionCreatedIn << "\n";
      outFile << "# Line : " << lineNumber << "\n";
      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      //outFile.precision(14);
      outFile.precision(std::numeric_limits<t>::digits10);
      for(unsigned row = 0; row < array.size(); ++row)
        {
          outFile << row << "\t";

          for(unsigned col = 0; col < array[row].size(); ++col)
            {
              if( col ==(array[row].size()-1) )
                {
                  outFile << array[row][col] << "\n";
                }
              else
                {
                  outFile << array[row][col] << "\t";
                }
            }

          // double space
          if((row+1)%togetherLines==0)
            {
              if(row !=(array.size()-1)){outFile << "\n\n";}
            }
        }
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      std::cerr << "could not open file " << fileName << std::endl;
      exit(EXIT_FAILURE);
    }

  outFile.close();

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t>
void write(const double* const array, const unsigned& numRows, const unsigned& numCols,
           const std::vector<std::string>& info, const unsigned& togetherLines)
{
  // the following will write the array but introduce double lines after
  // writing togetherLines number or lines if togetherLines=4 then there
  // will be a double space after every 4 lines, this is needed to use
  // the `using index' feature of gnuplot

  const std::string fileName = info[0];
  const std::string fileCreatedIn = info[1];
  const std::string functionCreatedIn = info[2];
  const std::string lineNumber = info[3];

  assert(!fileExists(fileName));

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good())
    {
      outFile << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      outFile << "# File : " << fileCreatedIn << "\n";
      outFile << "# Function : " << functionCreatedIn << "\n";
      outFile << "# Line : " << lineNumber << "\n";
      outFile << "# numRows x numCols :: " << numRows << " x " << numCols << "\n";
      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      for(unsigned row = 0; row < numRows; ++row)
        {
          outFile << row << "\t";

          for(unsigned col = 0; col < numCols; ++col)
            {
              if(col ==(numCols-1))
                {
                  outFile << array[row*numCols+col] << "\n";
                }
              else
                {
                  outFile << array[row*numCols+col] << "\t";
                }
            }

          // double space
          if((row+1)%togetherLines==0)
            {
              if(row ==(numRows-1))
                {
                }
              else
                {
                  outFile << "\n\n";
                }
            }
        }
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }

  outFile.close();

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t>
void write(t const value, const std::string& fileName){
  // write a single value to a file

  assert(!fileExists(fileName.c_str()));

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good()){
    outFile << "#" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
    outFile.precision(std::numeric_limits<t>::digits10);
    outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    outFile << value << "\n";
    outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  } else {
    std::cout << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
    std::cerr << alert("ERROR: ")
              << "could not write to stream " << std::endl;
    exit(EXIT_FAILURE);
  }

  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
template <typename t>
void write(const std::vector<t>& array, const std::string& fileName)
{
  // quick and dirty testing
  assert(!fileExists(fileName.c_str()));

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  const unsigned numElements = array.size();

  if(outFile.good())
    {
      outFile << "#" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      outFile << "# " << numElements << " elements " << std::endl;
      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      outFile.precision(std::numeric_limits<t>::digits10);
      for(unsigned col = 0; col < numElements; ++col)
        {
          if(col == (numElements-1))
            {
              outFile << array[col] << std::endl;
            }
          else
            {
              outFile << array[col] << "\t";
            }
        }
      outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR: ")
                << "could not write to steam " << std::endl;
      exit(EXIT_FAILURE);
    }

  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
template <class t>
void write(const t* array, const size_t size, const size_t stride, const std::string& fileName)
{
  assert(!fileExists(fileName.c_str()));
  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good())
    {
      outFile << "#"<< __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      outFile.precision(std::numeric_limits<t>::digits10);
      outFile << "# " << size << " elements " << std::endl;
      for(size_t i = 0; i < size; i += stride)
        {
          if(i == (size-1))
            {
              outFile << array[i] << std::endl;
            }
          else
            {
              outFile << array[i] << "\t";
            }
        }
      outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cout << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR: ") << "could not open file " << fileName << std::endl;
      exit(EXIT_FAILURE);
    }

  outFile.close();

  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
template <class t> void
write(const t* const array, const size_t size, const size_t stride, const size_t numRows,
      const size_t numCols, const std::string& fileName){
  // 2-dim array
  // this is used for quick and dirty testing

  assert(!fileExists(fileName.c_str()));

  assert(size >= numRows*numCols);

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good()){
    outFile << "#"<< __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
    outFile << "# " << numRows << " x " << numCols << "\n";

    outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    outFile.precision(std::numeric_limits<t>::digits10);

    for(size_t row = 0; row < numRows; ++row){
      outFile << row << "\t";
      for(size_t col = 0; col < numCols; ++col){
        if(col == (numCols-1)){
          outFile << array[stride*(row*numCols+col)] << std::endl;
        } else {
          outFile << array[stride*(row*numCols+col)] << "\t";
        }
      }
    }
    outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  } else {
    std::cout << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
    std::cerr << alert("ERROR: ") << "could not open file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  return;
}
//------------------------------




//------------------------------------------------------------------------------------------
template <typename t>
void write(const std::vector<std::vector<std::vector<t> > >& array,
           const std::vector<std::string>& info)
{
  const std::string fileName = info[0];
  const std::string fileCreatedIn = info[1];
  const std::string functionCreatedIn = info[2];
  const std::string lineNumber = info[3];

  assert(!fileExists(fileName));

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good())
    {
      outFile << "# " << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      outFile << "# File : " << fileCreatedIn << "\n";
      outFile << "# Function : " << functionCreatedIn << "\n";
      outFile << "# Line : " << lineNumber << "\n\n";
      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      //outFile.preion(14);
      outFile.precision(std::numeric_limits<t>::digits10);
      for(unsigned n = 0, table = 0; table < array.size(); ++table)
        {
          outFile << "# table " << table << "\n";
          if(!array[table].size()){outFile << "?\n" << "\n";}

          for(unsigned row = 0; row < array[table].size(); ++row)
            {
              outFile << n++ << "\t";
              for(unsigned col = 0; col < array[table][row].size(); ++col)
                {
                  if(col == (array[table][row].size()-1))
                    {
                      outFile << array[table][row][col] << "\n";
                    }
                  else
                    {
                      outFile << array[table][row][col] << "\t";
                    }
                }
              outFile << "\n";
            }

          if(table != (array.size()-1)){outFile << "\n";}

        }
      outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      std::cerr << alert("ERROR: ")
                << "could not open file " << fileName << std::endl;
      exit(EXIT_FAILURE);
    }

  outFile.close();

  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
template <typename t>
void write(const std::vector<std::vector<std::vector<t> > >& array, std::ofstream& ofs,
           const std::vector<std::string>& info)
{
  const std::string message = info[0];
  const std::string fileCreatedIn = info[1];
  const std::string functionCreatedIn = info[2];
  const std::string lineNumber = info[3];

  if(ofs.good())
    {
      ofs << "# " << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      ofs << "# File : " << fileCreatedIn << "\n";
      ofs << "# Function : " << functionCreatedIn << "\n";
      ofs << "# Line : " << lineNumber << "\n";
      ofs << "# " << message << "\n\n";
      ofs.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      //ofs.precision(14);
      ofs.precision(std::numeric_limits<t>::digits10);
      for(unsigned n = 0, table = 0; table < array.size(); ++table)
        {
          ofs << "# table " << table << "\n";
          if(!array[table].size()){ofs << "?\n" << "\n";}

          for(unsigned row = 0; row < array[table].size(); ++row)
            {
              ofs << n++ << "\t";
              for(unsigned col = 0; col < array[table][row].size(); ++col)
                {
                  if(col == (array[table][row].size()-1))
                    {
                      ofs << array[table][row][col] << "\n";
                    }
                  else
                    {
                      ofs << array[table][row][col] << "\t";
                    }
                }
              ofs << "\n";
            }

          //if(table != (array.size()-1)){ofs << "\n";}

        }
      ofs.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      std::cerr << alert("ERROR: ") << "could not open stream " << std::endl;
      exit(EXIT_FAILURE);
    }

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <typename t>
void write(const std::vector<std::vector<std::vector<t> > >& array,
           const std::string& fileName)
{
  //Mon Dec 20 23:17:37 2010 quick and dirty testing

  // check if the file exists
  assert(!fileExists(fileName));

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good())
    {
      outFile << "# " << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      //outFile.precision(14);
      outFile.precision(std::numeric_limits<t>::digits10);
      for(unsigned table = 0; table < array.size(); ++table)
        {
          outFile << "# table " << table << "\n";
          if(!array[table].size()){outFile << "?\n" << "\n";}

          for(unsigned row = 0; row < array[table].size(); ++row)
            {
              outFile << row << "\t";
              for(unsigned col = 0; col < array[table][row].size(); ++col)
                {
                  if(col == (array[table][row].size()-1))
                    {
                      outFile << array[table][row][col] << "\n";
                    }
                  else
                    {
                      outFile << array[table][row][col] << "\t";
                    }
                }
              outFile << "\n";
            }

          if(table !=(array.size()-1)){outFile << "\n";}
        }
      outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      std::cerr << alert("ERROR: ")
                << "could not open file " << fileName << std::endl;
      exit(EXIT_FAILURE);
    }

  outFile.close();

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <class t>
void write(const std::vector<t>& array, std::ofstream& ofs)
{
  if(ofs.good())
    {
      ofs.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      for(unsigned col = 0; col < array.size(); ++col)
        {
          if(col == (array.size()-1))
            {
              ofs << array[col] << std::endl;
            }
          else
            {
              ofs << array[col] << "\t";
            }
        }
      ofs.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR: ")
                << "could not write to steam " << std::endl;
      exit(EXIT_FAILURE);
    }

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <class t>
void write(const std::vector<std::vector<t> >& array, const std::string& fileName)
{
  // 2-dim
  // for quick and dirty testing during debugging phase
  std::cout << "\n" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";

  assert(!fileExists(fileName.c_str()));

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good())
    {
      outFile << "#"<< __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      outFile.precision(std::numeric_limits<t>::digits10);
      for(unsigned row = 0; row < array.size(); ++row)
        {
          if(array[row].size())
            {
              outFile << row << "\t";
            }
          else
            {
              outFile << row << "\n";
            }

          for(unsigned col = 0; col < array[row].size(); ++col)
            {
              if(col == (array[row].size()-1))
                {
                  outFile << array[row][col] << "\n";
                }
              else
                {
                  outFile << array[row][col] << "\t";
                }
            }
        }
      outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cerr << "File: " << __FILE__ << std::endl;
      std::cerr << "Function: " << __PRETTY_FUNCTION__ << std::endl;
      std::cerr << "ERROR: "
                << "could not open file " << fileName << std::endl;
      exit(EXIT_FAILURE);
    }

  outFile.close();

  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
template <class t>
void write(const t* const data, const unsigned& numRows, const unsigned& numCols,
           const std::vector<std::string>& info)
{
  const std::string fileName = info[0];
  const std::string fileCreatedIn = info[1];
  const std::string functionCreatedIn = info[2];
  const std::string lineNumber = info[3];

  if(fileExists(fileName))
    {
      std::cout << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR:: ")
                << fileName << " already exists ...\n";
      exit(EXIT_FAILURE);
    }

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good())
    {
      outFile << "# written by :" << __FILE__ << " --> " << __func__
              << "() --> " << __LINE__ << std::endl;
      outFile << "# File : " << fileCreatedIn << std::endl;
      outFile << "# Function : " << functionCreatedIn << std::endl;
      outFile << "# Line Number : " << lineNumber << std::endl;

      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      for(unsigned row = 0; row < numRows; ++row)
        {
          for(unsigned col = 0; col < numCols; ++col)
            {
              if(col == (numCols-1))
                {
                  outFile << data[row*numCols+col] << std::endl;
                }
              else
                {
                  outFile << data[row*numCols+col] << "\t";
                }
            }
        }
      outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR::") << "stream not good " << fileName << std::endl;
      exit(EXIT_FAILURE);
    }


  outFile.close();
  return;
}
//------------------------------

//------------------------------------------------------------------------------------------
template <class t>
void write(const std::vector<std::vector<t> >& array, const std::vector<std::string>& info)
{
  const std::string fileName = info[0];
  const std::string fileCreatedIn = info[1];
  const std::string functionCreatedIn = info[2];
  const std::string lineNumber = info[3];

  if(fileExists(fileName))
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR:: ")
                << fileName << " already exists ...\n";
      exit(EXIT_FAILURE);
    }

  std::ofstream outFile(fileName.c_str(), std::ios::out);

  if(outFile.good())
    {
      outFile << "#"<< __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      outFile << "# File : " << fileCreatedIn << std::endl;
      outFile << "# Function : " << functionCreatedIn << std::endl;
      outFile << "# Line : " << lineNumber << std::endl;
      outFile.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      //outFile.precision(14);
      outFile.precision(std::numeric_limits<t>::digits10);
      for(unsigned row = 0; row < array.size(); ++row)
        {
          if(array[row].size())
            {
              outFile << row << "\t";
            }
          else
            {
              outFile << row << "\n";
            }

          for(unsigned col = 0; col < array[row].size(); ++col)
            {
              if(col == (array[row].size()-1))
                {
                  outFile << array[row][col] << std::endl;
                }
              else
                {
                  outFile << array[row][col] << "\t";
                }
            }
        }
      outFile.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cout << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR: ")
                << "could not open file " << fileName << std::endl;
      exit(EXIT_FAILURE);
    }

  outFile.close();

  return;
}
//------------------------------


//------------------------------------------------------------------------------------------
template <class t>
void write(const std::vector<std::vector<t> >& array, std::ofstream& ofs,
           const std::vector<std::string>& info)
{
  const std::string message = info[0];
  const std::string fileCreatedIn = info[1];
  const std::string functionCreatedIn = info[2];
  const std::string lineNumber = info[3];

  if(ofs.good())
    {
      ofs << "#"<< __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      ofs << "# File : " << fileCreatedIn << "\n";
      ofs << "# Function : " << functionCreatedIn << "\n";
      ofs << "# Line : " << lineNumber << "\n";
      ofs << "# " << message << "\n\n";
      ofs.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
      //ofs.precision(14);
      ofs.precision(std::numeric_limits<t>::digits10);
      for(unsigned row = 0; row < array.size(); ++row)
        {
          if(array[row].size())
            {
              ofs << row << "\t";
            }
          else
            {
              ofs << row << "\n";
            }

          for(unsigned col = 0; col < array[row].size(); ++col)
            {
              if(col == (array[row].size()-1))
                {
                  ofs << array[row][col] << "\n";
                }
              else
                {
                  ofs << array[row][col] << "\t";
                }
            }
        }
      ofs.unsetf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
    }
  else
    {
      std::cerr << __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
      std::cerr << alert("ERROR: ")
                << "could not open stream " << std::endl;
      exit(EXIT_FAILURE);
    }

  return;
}
//------------------------------
#endif
