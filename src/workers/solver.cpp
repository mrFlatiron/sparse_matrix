#include "sparse_matrix/msr_pthread.h"
#include "containers/cycle_buf.h"
#include <pthread.h>
#include <cstdio>

void *solve (void *args)
{
  msr_pthread_initializer *initer = (msr_pthread_initializer *)args;
  msr_pthread handler (*initer);

  printf ("thread - %d\n", handler.t_id ());
  std::vector<double> in = {1, 2, 3, 4, 5};
  std::vector<double> out (5);

  handler.mult_vector (in, out);

  if (handler.t_id () == 0)
    {
      for (int i = 0; i < 5; i++)
        printf ("out[%d] = %.3lf\n", i, out[i]);
    }


}
