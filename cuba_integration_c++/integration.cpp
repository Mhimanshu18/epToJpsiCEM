#include <iostream>
#include <cmath>

#include "cuba.h"

static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{

  ff[0] = xx[0];

    return 0;
}

/*********************************************************************/

#define NDIM 1
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 2
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 50000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

int main()
{
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
        EPSREL, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        GRIDNO, STATEFILE, SPIN,
        &neval, &fail, integral, error, prob);

  std::cout << (double)integral[NCOMP] << "\t" << (double)error[NCOMP] << "\t" << (double)prob[NCOMP] << "\n";

  return 0;
}
