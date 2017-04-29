#include "sparse_matrix/msr_thread_dqgmres_solver.h"
#include "containers/cycle_buf.h"
#include "threads/thread_vector_utils.h"
#include "sparse_matrix/msr_matrix.h"
#include <pthread.h>
#include <cstdio>

void *test (void *args)
{
  msr_thread_dqgmres_solver handler (*((msr_thread_dqgmres_solver *)args));

  int n = handler.n ();

  printf ("thread - %d\n", handler.t_id ());
  std::vector<double> in = {1, 2, 3, 4, 5};
  std::vector<double> out (n);

  handler.mult_vector (handler.matrix (), in, out);

  if (handler.t_id () == 0)
    {
      for (int i = 0; i < n; i++)
        printf ("out[%d] = %.3lf\n", i, out[i]);
      printf ("out - 2 * in = \n");
    }
  std::vector<double> res (n);
  thread_utils::lin_combination_1 (handler, handler.shared_ref (), out, res, in, -2);
  if (handler.t_id () == 0)
    {
      for (int i = 0; i < n; i++)
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

  thread_utils::l2_norm (handler, in, handler.p_sized_buf ());

  if (handler.t_id () == 0)
      printf ("l2")

    return args;
}

void *solve(void *args)
{
    msr_thread_dqgmres_solver handler (*((msr_thread_dqgmres_solver *)args));

    handler.dqgmres_solve ();
    return args;
}
