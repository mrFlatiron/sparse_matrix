#include "sparse_matrix/msr_matrix.h"
#include "sparse_matrix/msr_pthread.h"
#include <pthread.h>
#include "workers/solver.h"
#include "containers/cycle_buf.h"
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
  std::vector<double> buf (5);

  msr_pthread_initializer *initers = new msr_pthread_initializer[p];

  for (int t = 1; t < p; t++)
    {

      initers[t].matrix = &msr;
      initers[t].p = p;
      initers[t].buf = &buf;
      initers[t].barrier = &barrier;
      initers[t].t = t;

      pthread_create (&pt, NULL, solve, initers + t);
    }
  initers[0].matrix = &msr;
  initers[0].p = p;
  initers[0].buf = &buf;
  initers[0].barrier = &barrier;
  initers[0].t = 0;
  solve (initers + 0);

  cycle_buf<int> test_buf (5);
  for (int i = 1; i < 9; i++)
    {
      test_buf.insert (i);

      test_buf.to_start ();

      int k = 0;
      while (!test_buf.is_loop_done ())
        {
          int j = test_buf.get_next ();
          printf ("buf[%d] = %d\n", k, j);
          k++;
        }
      printf ("------------\n");
    }
  return 0;
}
