#include <iostream>
#include <fstream>
#include <cstdio>

#include <cstdlib> // for atof()

//#define debug

#include <cstring> // for strtok()

#include <vector>

#include <string>

#include "ap_utils.h"

// void read(std::vector<std::vector<double> > & data, std::string const filename) {
//   std::cout << "\n" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
//   std::ifstream ifs(filename.c_str(), std::ios::in);
//   std::string line;
//   while(std::getline(ifs, line)) {
  

//------------------------------------------------------------------------------------------
//Sat Aug 28 01:35:31 2010
void read(const std::string & fileName, std::vector <std::vector <double> > & data, 
	  const unsigned & numRows, const unsigned & numCols){
  // read an array from a file
  
  // delimters are hash, space, comma, tab
  char delimters[] = " ,\t" ;

  double tempDbl = 0.0;

  char* tempChr = NULL;

  std::ifstream file(fileName.c_str());
  if(!file.good())
    {
      std::cerr << __FILE__ << " --> " << __FUNCTION__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR");
      exit(EXIT_FAILURE);
    }

  // assuming each line of data file 
  // does not exceed 512 characters
  // an extra character for the null character that is automatically
  // appended 
  char line[512+1];

  unsigned row = 0;
  unsigned col = 0;

  while(!file.getline(line, sizeof(line), '\n').eof()) {
      col = 0 ;
      tempChr = strtok(line, delimters);

      if(!tempChr){//if tempChr == NULL
	  continue;
	}

      if(tempChr[0] == '\n'){ continue; }

      if(tempChr[0] == '#'){ continue; }

      tempDbl = atof(tempChr);

      data[row][col] = tempDbl;

      for(col = 1; col < numCols; ++col){
	  tempDbl = atof(strtok(NULL, delimters));
	  data[row][col] = tempDbl;
	}

      row = row + 1;
      SHOW(row);
      if(row == (numRows+1)){ break; }
    }

  // close the input file
  file.close();

  return ;
}
//------------------------------



//------------------------------------------------------------------------------------------
void read(const std::string & fileName, double* const data, 
	  const unsigned numRows, const unsigned numCols) 
{
  // delimters are hash, space, comma, tab
  char delimters[] = " ,\t";

  double tempDbl = 0.0;

  char* tempChr = nullptr;

  std::FILE* file = nullptr; 

  file = fopen(fileName.c_str(), "r");
  if (!file) {
      std::cerr << __FILE__ << " --> " << __FUNCTION__ << "() --> " << __LINE__ << std::endl;
      std::cerr << alert("ERROR :: ") << "could not open " << fileName << std::endl; 
      exit(EXIT_FAILURE);
    }

  // assuming each line of data file 
  // does not exceed 512 characters
  // an extra character for the null character that is automatically
  // appended 
  char line[512+1];

  unsigned int row = 0;
  unsigned int col = 0;

  while( fgets(line, sizeof(line), file) )    {
    //std::cout << line << std::endl;
      col = 0;
      tempChr = strtok(line, delimters);

      if (tempChr[0] == '\n' || tempChr[0] == '#') { continue; }

      tempDbl = atof(tempChr);

      data[row*numCols+col] = tempDbl;

      for(col = 1; col < numCols; ++col)
	{
	  tempDbl = atof(strtok(NULL, delimters));
	  data[row*numCols+col] = tempDbl;
	}

      row = row + 1;
      if( row == numRows )
	{
	  break;
	}
    }

  fclose(file);

  return;
}
