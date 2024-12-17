#include <iostream>
#include <cmath>

#include "cuba.h"

double myFunction(double x, double y)
{

  return exp(x) * sin(y);
}

static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
  double xmin = 2.0, xmax = 5.0;
  double dx = (xmax - xmin);
  double x = xmin + dx * xx[0];

  double ymin = 1.25, ymax = 3.0;
  double dy = (ymax - ymin);
  double y = ymin + dy * xx[1];

  ff[0] = dx * dy * myFunction(x, y);

  return 0;
}

#define NDIM 2
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 500000
#define MAXEVAL 5000000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

int main()
{
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
        EPSREL, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        GRIDNO, STATEFILE, SPIN,
        &neval, &fail, integral, error, prob);

  std::cout << (double)integral[0] << "\t" << (double)error[0] << "\t" << (double)prob[0] << "\t" << fail << "\n";

  return 0;
}
