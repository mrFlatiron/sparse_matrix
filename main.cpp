#include "sparse_matrix/msr_matrix.h"
#include "sparse_matrix/msr_pthread.h"
#include "pthread.h"
#include "workers/solver.h"
#include <cstdlib>

using namespace std;

int main (int argc, char *argv[])
{
  if (argc < 2)
    {
      printf ("Usage : %s p\n", argv[0]);
      return 0;
    }
  int p = atoi (argv[1]);
  if (p <= 0)
    {
      printf ("Wrong p\n");
      return 0;
    }
  msr_matrix msr;
  msr.convert (5, {1, 0, 0, 2, 0,
                   3, 4, 0, 0, 0,
                   0, 0, 7, 0, 0,
                   8, 9, 10, 11, 0,
                   0, 0, 0, 0, 12});
  msr.dump ();

  pthread_barrier_t barrier;
  pthread_barrier_init (&barrier, NULL, p);
  pthread_t pt;
  double *buf = new double[5];

  msr_pthread_initializer *initers = new msr_pthread_initializer[p];

  for (int t = 1; t < p; t++)
    {

      initers[t].matrix = &msr;
      initers[t].p = p;
      initers[t].buf = buf;
      initers[t].barrier = &barrier;
      initers[t].t = t;

      pthread_create (&pt, NULL, solve, initers + t);
    }
  initers[0].matrix = &msr;
  initers[0].p = p;
  initers[0].buf = buf;
  initers[0].barrier = &barrier;
  initers[0].t = 0;
  solve (initers + 0);
  return 0;
}
