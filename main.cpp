#include "sparse_matrix/msr_matrix.h"
#include "sparse_matrix/msr_thread_dqgmres_solver.h"
#include <pthread.h>
#include "workers/solver.h"
#include "containers/cycle_buf.h"
#include <cstdlib>

using namespace std;



int main (int argc, char *argv[])
{
  const int max_iter = 100;
  const double stop_criterion = 1e-6;
  const int n = 5;

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
      printf ("Wrong dim\n");
      return 0;
    }
  msr_matrix msr;
  msr.convert (n, {1, 0, 0, 0, 0,
                   0, 2, 3, 0, 0,
                   0, 4, 1, 0, 0,
                   0, 2, 1, 3, 0,
                   0, 0, 0, 2, 1});

  msr_matrix save;
  save = msr;
  msr.dump ();

  std::vector<double> x (n);

  for (int i = 0; i < n; i++)
    x[i] = i & 1;

  msr_matrix precond;

  pthread_barrier_t barrier;
  pthread_barrier_init (&barrier, NULL, p);
  pthread_t pt;
  std::vector<double> buf (n);
  std::vector<double> v1 (n, 0), v2 (n), v3 (n), p_sized (p);
  std::vector<double> rhs (n);
  msr.mult_vector (x, rhs);
  std::vector<double> rhs_save (rhs);

  cycle_buf<std::vector<double>> basis (dim);
  cycle_buf<std::vector<double>> basis_derivs (dim);
  cycle_buf<std::vector<double>> turns (dim);
  std::vector<double> hessenberg (dim + 2);

  std::vector<msr_thread_dqgmres_solver> handlers;
  bool flag = false;
  for (int i = 0; i < p; i++)
    {
      handlers.push_back (msr_thread_dqgmres_solver
                          (i, p, &barrier, buf, msr, precond,
                           preconditioner_type::identity, dim,
                           max_iter, stop_criterion, flag, rhs,
                           basis, basis_derivs, turns,
                           hessenberg,
                           p_sized, v1, v2, v3));
    }


  for (int t = 1; t < p; t++)
    {
      pthread_create (&pt, NULL, solve, handlers.data () + t);
    }
  solve (handlers.data () + 0);

  std::vector<double> rhs_comp (n);
  save.mult_vector (v1, rhs_comp);

  for (int i = 0; i < n; i++)
    printf ("rhs_comp[%d] = %lf\n", i, rhs_comp[i] - rhs_save[i]);

  return 0;
}
