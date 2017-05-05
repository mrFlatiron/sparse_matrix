#include "sparse_matrix/msr_matrix.h"
#include "sparse_matrix/msr_thread_dqgmres_solver.h"
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
  const double stop_criterion = 1e-15;
  const int n = 4000;

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
//  msr.convert (n, {1, 0, 0, 0, 0,
//                   0, 2, 4, 2, 0,
//                   0, 4, 1, 1, 0,
//                   0, 2, 1, 3, 0,
//                   0, 0, 0, 2, 1});

//  msr.convert (n, {3, 0, 5, 0, 10,
//                   0, 2, 0, 4, 0,
//                   2, 4, 4, 0, 5,
//                   0, 0, 0, 5, 0,
//                   1, 0, 2, 0, 6});

//  msr.convert (n, {1, 0, 0, 0, 0,
//                   0, 2, 3, 0, 0,
//                   0, 0, 1, 2, 0,
//                   0, 0, 0, 3, 0,
//                   0, 0, 0, 0, 1});

//  msr.convert (n,     {1, 2, 0 , 0, 0, 3, 0, 0, 0, 0,
//                       2, 1, 2 , 0, 0, 0, 3, 0, 0, 0,
//                       0, 2, 1 , 2, 0, 0, 0, 3, 0, 0,
//                       0, 0, 2 , 1, 2, 0, 0, 0, 3, 0,
//                       0, 0, 0 , 2, 1, 2, 0, 0, 0, 3,
//                       3, 0, 0 , 0, 2, 1, 2, 0, 0, 0,
//                       0, 3, 0 , 0, 0, 2, 1, 2, 0, 0,
//                       0, 0, 3 , 0, 0, 0, 2, 1, 2, 0,
//                       0, 0, 0 , 3, 0, 0, 0, 2, 1, 2,
//                       0, 0, 0 , 0, 3, 0, 0, 0, 2, 1});


  simple_vector big_matrix (n * n);
  srand (1);
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
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

          if (j == i + n / 4)
            {
              big_matrix[i * n + j] = 2;
              continue;
            }
          if (j == i + n / 4  + 1)
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

  msr.convert (n, big_matrix);

  msr_matrix save;
  save = msr;
  msr.dump ();

  simple_vector x (n);

  for (int i = 0; i < n; i++)
    x[i] = 1;

  msr_matrix precond;

  pthread_barrier_t barrier;
  pthread_barrier_init (&barrier, NULL, p);
  pthread_t pt;
  simple_vector buf (n);
  simple_vector x_inout (n), *v1, v2 (n), p_sized (p);
  simple_vector rhs (n);
  msr.mult_vector (x, rhs);
  simple_vector rhs_save (rhs);
  simple_vector rhs_true_save (rhs);

  limited_deque<simple_vector> basis (dim);
  limited_deque<simple_vector> basis_derivs (dim);
  limited_deque<simple_vector> turns (dim);
  simple_vector hessenberg (dim + 2);

  std::vector<msr_thread_dqgmres_solver> handlers;
  bool flag = false;
  for (int i = 0; i < p; i++)
    {
      handlers.push_back (msr_thread_dqgmres_solver
                          (i, p, &barrier, buf, msr, precond,
                           preconditioner_type::jacobi, dim,
                           max_iter, stop_criterion, flag, rhs,
                           rhs_save,
                           basis, basis_derivs, turns,
                           hessenberg,
                           p_sized, x_inout, &v1, v2));
    }


  for (int t = 1; t < p; t++)
    {
      pthread_create (&pt, NULL, solve, handlers.data () + t);
    }
  solve (handlers.data () + 0);

  simple_vector rhs_comp (n);
  save.mult_vector (x_inout, rhs_comp);

  double res = 0;
  for (int i = 0; i < n; i++)
    {
      res += pow ((rhs_comp[i] - rhs_true_save[i]), 2);
    }
  res = sqrt (res);
  printf ("Residual L2 : %.6le\n", res);
  return 0;
}
