#include "sparse_matrix/msr_matrix.h"
#include "sparse_matrix/msr_thread_dqgmres_solver.h"
#include "sparse_matrix/msr_dqgmres_initializer.h"
#include <pthread.h>
#include "workers/solver.h"
#include "containers/limited_deque.h"
#include <cstdlib>
#include <cmath>
#include "containers/simple_vector.h"

using namespace std;



int main (int argc, char *argv[])
{
  const int max_iter = 1000;
  const double stop_criterion = 1e-9;
  const int n = 6000;

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

  simple_vector big_matrix (n * n);

  srand (1);
  for (int i = 0; i < n; i++)
    {
      for (int j = i; j < n; j++)
        {
          if (i == j)
            {
              big_matrix[i * n + j] = 7;
              continue;
            }
          if (j == i + 1)
            {
              big_matrix[i * n + j] = 1;
              continue;
            }

          if (j == i + n / 3)
            {
              big_matrix[i * n + j] = 2;
              continue;
            }
          if (j == i + n / 3  + 1)
            {
              big_matrix[i * n + j] = 3;
              continue;
            }
          big_matrix[i * n + j] = 0;
        }
    }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      big_matrix[i * n + j] = big_matrix[j * n + i];

//  for (int i = 0; i < n; i++)
//    for (int j = 0; j < n; j++)
//      big_matrix[i * n + j] = 1. / (1 + i + j);

  msr.convert (n, big_matrix);
  pthread_t pt;
  if (n <= 20)
  msr.dump ();

  simple_vector x (n);

  for (int i = 0; i < n; i++)
    x[i] = 10;

  simple_vector rhs (n);
  msr.mult_vector (x, rhs);

  simple_vector x_ini (n);
  for (int i = 0; i < n ; i++)
    x_ini[i] = 0;

  simple_vector x_out (n);

  msr_dqgmres_initializer initializer (p, msr,
                                       preconditioner_type::jacobi,
                                       dim, max_iter, stop_criterion,
                                       x_ini, x_out, rhs);

  for (int i = 0; i < n; i++)
    initializer.m_v2[i] = 2;

  std::vector<msr_thread_dqgmres_solver> handlers;
  for (int i = 0; i < p; i++)
    {
      handlers.push_back (msr_thread_dqgmres_solver (i, initializer));
    }


  for (int t = 1; t < p; t++)
    {
      pthread_create (&pt, NULL, solve, handlers.data () + t);
    }
  solve (handlers.data () + 0);

  simple_vector rhs_comp (n);
  msr.mult_vector (x_out, rhs_comp);

  double res = 0;
  double desc = 0;
  for (int i = 0; i < n; i++)
    {
      desc += pow (rhs_comp[i] - rhs[i], 2);
      res += pow (x_out[i] - x[i], 2);
    }
  res = sqrt (res);
  printf ("Residual L2 : %.6le, Discrepancy L2 : %.6le\n", res, desc);
  return 0;
}
