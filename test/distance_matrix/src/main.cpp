#include <iostream>
#include <fstream>
#include <vector>

#include <algorithm>

#include <cstdlib>

#include "ap_utils.h"


//------------------------------------------------------------------------------------------
int test() {
  srand(time(0));
  
  const unsigned numData = 1;
  const unsigned numCenters = 100;
  const unsigned ndims = 30;
  
  std::vector<std::vector<double> > dm(numData, std::vector<double>(numCenters));

  std::vector<std::vector<double> > data(numData, std::vector<double>(ndims));
  std::for_each(data.begin(), data.end(), [](std::vector<double>& v)
                {
                  std::generate(v.begin(), v.end(), [](){return rand()%10;});
                });
  //display(data);
  
  std::vector<std::vector<double> > centers(numCenters, std::vector<double>(ndims));
  std::for_each(centers.begin(), centers.end(), [](std::vector<double>& v)
                {
                  std::generate(v.begin(), v.end(), [](){return rand()%10;});
                });
  //display(centers);

  distanceMatrix(dm, data, centers);

  //display(dm);

  write(data, "d.txt");
  write(centers, "c.txt");
  write(dm, "dc.txt");

  return 0 ;
}
//------------------------------


//------------------------------------------------------------------------------------------
int test1(size_t const numdata, size_t const numcenters, size_t const ndims)
{
  double* d = new double [numdata*ndims];

  double* c = new double [numcenters*ndims];
      
  {
    gsl_rng* grng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(grng, time(NULL));

    std::for_each(d, d+numdata*ndims, [&grng](double&i){i = gsl_rng_uniform(grng);});
    std::for_each(c, c+numcenters*ndims, [&grng](double&i){i = gsl_rng_uniform(grng);});

    gsl_rng_free(grng);
  }

  double* z = new double [numdata*numcenters];

  distanceMatrix(z, d, numdata, c, numcenters, ndims);

  {
    write(d, numdata*ndims, 1, numdata, ndims, "d.out");
    write(c, numcenters*ndims, 1, numcenters, ndims, "c.out");
    write(z, numdata*numcenters, 1, numdata, numcenters, "z.out");
    
    FILE* oct = popen("/usr/bin/octave -q", "w");
    fprintf(oct, "load d.out; load c.out; load z.out;\n");
    fprintf(oct, "d=d(:,2:end); c=c(:,2:end); z=z(:,2:end);\n");
    fprintf(oct, "dm = sqrt(squeeze(sumsq(d - permute(c, [3 2 1]), 2)));\n");
    fprintf(oct, "norm(z-dm)");
    fclose(oct);
  }

  delete [] d;
  delete [] c;
  delete [] z;
  
  return 0;
}
//------------------------------


//------------------------------------------------------------------------------------------
int main()
{
  test1(10, 10, 3);

  return 0;
}
//------------------------------
