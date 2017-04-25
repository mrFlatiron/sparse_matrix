#include "sparse_matrix/msr_matrix.h"
#include "sparse_matrix/msr_pthread.h"
#include <pthread.h>
#include "workers/solver.h"
#include "containers/cycle_buf.h"
#include <cstdlib>

using namespace std;

const int max_iter = 100;
const double stop_criterion = 1e-6;

int main (int argc, char *argv[])
{
  if (argc < 3)
    {
      printf ("Usage : %s p, dim\n", argv[0]);
      return 0;
    }
  int p = atoi (argv[1]);
  if (p <= 0)
    {
      printf ("Wrong p\n");
      return 0;
    }

  int dim = atoi (argv[2]);
  if (dim <= 0)
    {
      printf ("Wrong d\n");
      return 0;
    }
  msr_matrix msr;
  msr.convert (5, {1, 0, 0, 2, 0,
                   3, 4, 0, 0, 0,
                   0, 0, 7, 0, 0,
                   8, 9, 10, 11, 0,
                   0, 0, 0, 0, 12});
  msr.dump ();

  msr_matrix precond;

  pthread_barrier_t barrier;
  pthread_barrier_init (&barrier, NULL, p);
  pthread_t pt;
  std::vector<double> buf (5);
  std::vector<double> v1 (5), v2 (5), v3 (5), p_sized (p);
  std::vector<double> rhs ({1, 2, 3, 4, 5});

  cycle_buf<std::vector<double>> basis (5);
  cycle_buf<std::vector<double>> basis_derivs (5);
  cycle_buf<std::vector<double>> turns (5);

{
  std::vector<msr_thread_dqgmres_solver> handlers;

  for (int i = 0; i < p; i++)
    {
      handlers.push_back (msr_thread_dqgmres_solver
                          (i, p, &barrier, buf, msr, precond,
                           preconditioner_type::jacobi, dim,
                           max_iter, stop_criterion, rhs,
                           basis, basis_derivs, turns,
                           p_sized, v1, v2, v3));
    }


  for (int t = 1; t < p; t++)
    {
      pthread_create (&pt, NULL, solve, handlers.data () + t);
    }
  solve (handlers.data () + 0);
}
  return 0;
}
