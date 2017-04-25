#include "sparse_matrix/msr_pthread.h"
#include "containers/cycle_buf.h"
#include "threads/thread_vector_utils.h"
#include <pthread.h>
#include <cstdio>

void *solve (void *args)
{
  msr_thread_dqgmres_solver handler (*((msr_thread_dqgmres_solver *)args));

  printf ("thread - %d\n", handler.t_id ());
  std::vector<double> in = {1, 2, 3, 4, 5};
  std::vector<double> out (5);

  handler.mult_vector (handler.matrix (), in, out);

  if (handler.t_id () == 0)
    {
      for (int i = 0; i < 5; i++)
        printf ("out[%d] = %.3lf\n", i, out[i]);
      printf ("out - 2 * in = \n");
    }
  std::vector<double> res (5);
  thread_utils::lin_combination_1 (handler, handler.shared_ref (), out, res, in, -2);
  if (handler.t_id () == 0)
    {
      for (int i = 0; i < 5; i++)
        printf ("res[%d] = %.3lf\n", i, res[i]);
      printf ("preconditioner:\n");
    }

  handler.compute_preconditioner ();

  if (handler.t_id () == 0)
    {
      handler.precond ().dump ();
    }

  handler.apply_preconditioner ();

  if (handler.t_id () == 0)
    {
      handler.matrix ().dump ();
    }


}
