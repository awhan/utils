#include <iostream>
#include <vector>

#include <cfloat>
// for dbl_epsilon

#include <gsl/gsl_math.h>
// for gsl_fcmp()

#include "ap_utils.h"

// #define debugCouts

// void dominanceCompare(const vector <double> & A, const vector <double> & B, 
// 		      const unsigned & size)

// A dominates B :: return 1
// B dominates A :: return -1
// A = B or A and B are non-dominant :: return 0
int dominanceCompare(const double * const A, const double * const B,
		     const unsigned & size)
{
  int result = 0;

  // not worse than, in all objectives
  bool conditionOne = false;

  // condition 1
  // A_i <= B_i for all i
  // if condition 1 is false then A worse
  // than B in one or more objective


  // strictly better than, in atleast one objective
  bool conditionTwo = false;

  // condition 2
  // there exists atleast one j for which
  // A_j < B_j 
  // that is A is strictly better than B in atleast
  // one objective



  // assume condition1 is true
  conditionOne = true; 
  for(unsigned i = 0; i < size; ++i)
    {
      //      if( A[i] > B[i] )
      if( gsl_fcmp(A[i], B[i], DBL_EPSILON) == 1 )
	{
	  conditionOne = false;
	  break;
	}
    }

	      
  // assume conditionTwo is not satisfied
  conditionTwo = false; 
  for(unsigned n = 0; n < size; ++n)
    {
      //      if( A[n] < B[n] )
      if( gsl_fcmp(A[n], B[n], DBL_EPSILON) == -1 )
	{
	  conditionTwo = true;
	  break;
	}
    }
  
#ifdef debug

  if ( conditionOne == true )
    {
      cout << "Condition 1 is true " << endl;
    }
  else
    {
      cout << "Condition 1 is false " << endl;
    }


  if ( conditionTwo == true )
    {
      cout << "Condition 2 is true " << endl;
    }
  else
    {
      cout << "Condition 2 is false " << endl;
    }
#endif


  // A dominates B :: return 1
  if( conditionOne == true && conditionTwo == true )
    {
      result = 1;

#ifdef debug
      cout << "(";
      for(unsigned n = 0; n < size-1; ++n)
	{
	  cout << A[n] <<", ";
	}
      cout << A[size-1] << ") " << endl 
	   << " dominates " << endl;

      cout << "(";
      for(unsigned n = 0; n < size-1; ++n)
	{
	  cout << B[n] <<", ";
	}
      cout << B[size-1] << ") " << endl;
#endif
    }
  // B dominates A :: return -1
  else if( conditionOne==false && conditionTwo==false )
    {
      result = -1;

#ifdef debug
      cout << "(";
      for(unsigned n = 0; n < size-1; ++n)
	{
	  cout << B[n] <<", ";
	}
      cout << B[size-1] << ") " << endl  
	   << " dominates " << endl;

      cout << "(";
      for(unsigned n = 0; n < size-1; ++n)
	{
	  cout << A[n] <<", ";
	}
      cout << A[size-1] << ") " << endl;
#endif
    }
  else
    // A = B or A and B are non-dominant
    {
      result = 0;
#ifdef debug
      cout << "either same points OR a non-dominant pair !!!" << endl;
#endif
    }

  return result;
}
//------------------------------


//------------------------------------------------------------------------------------------
int dominanceCompare(const std::vector <double> & a, const std::vector <double> & b) 
{
  // A dominates B :: return 1
  // B dominates A :: return -1
  // A = B or A and B are non-dominant :: return 0
  
  assert(a.size() == b.size());

  int result = 0;

  const unsigned size = a.size();


  // not worse than, in all objectives
  bool conditionOne = false;

  // condition 1
  // A_i <= B_i for all i
  // if condition 1 is false then A worse
  // than B in one or more objective


  // strictly better than, in atleast one objective
  bool conditionTwo = false;

  // condition 2
  // there exists atleast one j for which
  // A_j < B_j 
  // that is A is strictly better than B in atleast
  // one objective


  // assume condition1 is true
  conditionOne = true; 
  for(unsigned i = 0; i < size; ++i)
    {
      //      if( a[i] > b[i] )
      if(gsl_fcmp(a[i], b[i], DBL_EPSILON) == 1)
	{
	  conditionOne = false;
	  break;
	}
    }

	      
  // assume conditionTwo is not satisfied
  conditionTwo = false; 
  for(unsigned n = 0; n < size; ++n)
    {
      //      if( a[n] < b[n] )
      if(gsl_fcmp(a[n], b[n], DBL_EPSILON) == -1)
	{
	  conditionTwo = true;
	  break;
	}
    }
  
#ifdef debug

  if(conditionOne)
    {
      cout << "Condition 1 is true " << endl;
    }
  else
    {
      cout << "Condition 1 is false " << endl;
    }


  if(conditionTwo)
    {
      cout << "Condition 2 is true " << endl;
    }
  else
    {
      cout << "Condition 2 is false " << endl;
    }
#endif


  // A dominates B :: return 1
  if(conditionOne && conditionTwo)
    {
      result = 1;

#ifdef debug
      cout << "(";
      for(unsigned n = 0; n < size-1; ++n)
	{
	  cout << a[n] <<", ";
	}
      cout << a[size-1] << ") " << endl 
	   << " dominates " << endl;

      cout << "(";
      for(unsigned n = 0; n < size-1; ++n)
	{
	  cout << b[n] <<", ";
	}
      cout << b[size-1] << ") " << endl;
#endif
    }
  // B dominates A :: return -1
  else if(!conditionOne && !conditionTwo)
    {
      result = -1;

#ifdef debug
      cout << "(";
      for(unsigned n = 0; n < size-1; ++n)
	{
	  cout << b[n] <<", ";
	}
      cout << b[size-1] << ") " << endl  
	   << " dominates " << endl;

      cout << "(";
      for(unsigned n = 0; n < size-1; ++n)
	{
	  cout << a[n] <<", ";
	}
      cout << a[size-1] << ") " << endl;
#endif
    }
  else
    // A = B or A and B are non-dominant
    {
      result = 0;
#ifdef debug
      cout << "either same points OR a non-dominant pair !!!" << endl;
#endif
    }

  return result;
}


//------------------------------------------------------------------------------------------
bool nonDominated(const std::vector<double>& v, const std::vector<std::vector<double> >& a) 
{
  // is f non-dominated?
  for(auto itr = a.begin(); itr != a.end(); ++itr)
    {
      if(dominanceCompare(v, *itr)) return false;
    }
  return true;
}
//------------------------------


//------------------------------------------------------------------------------------------
void paretoSort(const std::vector<std::vector<double> >& data,
		std::vector<unsigned>& best)
{
  // takes the candidates for gapFillers and finds out the best
  // non-dominated set amongst them
  
  const unsigned numVec = data.size();
  
  if(numVec != best.size())
    {
      std::cerr << "File: " << __FILE__ << std::endl;
      std::cerr << "Function: " << __FUNCTION__  << std::endl;
      std::cerr << "Line: " << __LINE__  << std::endl;
      std::cerr << alert("ERROR:: ") 
		<< "size mismatch " << std::endl;   
      exit(EXIT_FAILURE);
    }
  
  for(unsigned n = 0; n < numVec; ++n)
    {
      best[n] = n;
    }
  
  std::vector<unsigned>::iterator iter;
  
  for(iter = best.begin(); iter != best.end(); ++iter)
    {
      const unsigned one = *iter;

      for(unsigned two = (one+1); two < numVec; ++two)
	{
	  const int result = dominanceCompare(data[one], data[two]);
	  
	  if(result==-1)
	    {
	      // two dominates one
	      // remove one from best
	      best.erase(iter);
	      //erase invalidates all iterator and references to elements after iter.
	      iter--;
	      break;
	    }
	  else if(!result)
	    {
	      // both non-dominated
	    }
	  else if(result==1)
	    {
	      // one dominates two
	      // remove two from best
	      best.erase(std::remove(best.begin(), best.end(), two), best.end());
	    }
	  else
	    {
	      std::cerr << "this should not happen" << std::endl;
	      break;
	    }
	}
    }

  return;
}
//------------------------------
