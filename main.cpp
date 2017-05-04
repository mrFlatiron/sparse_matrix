#include "sparse_matrix/msr_matrix.h"
#include "sparse_matrix/msr_thread_dqgmres_solver.h"
#include <pthread.h>
#include "workers/solver.h"
#include "containers/limited_deque.h"
#include <cstdlib>
#include <cmath>

using namespace std;



int main (int argc, char *argv[])
{
  const int max_iter = 10000;
  const double stop_criterion = 1e-3;
  const int n = 100;

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

//  msr.convert (n,     {1, 2, 3 , 0, 0, 0, 4, 5, 6, 0,
//                       0, 2, 3 , 0, 0, 1, 5, 4, 0, 0,
//                       0, 0, 2 , 1, 1, 0, 0, 5, 0, 3,
//                       2, 0, -1, 1, 2, 0, 0, 3, 1, 0,
//                       1, 0, 0 , 0, 1, 0, 0, 0, 0, 0,
//                       0, 2, 0 , 2, 0, 1, 0, 0, 2, 0,
//                       0, 0, 2 , 1, 0, 0, 1, 0, 0, 1,
//                       0, 2, 0 , 0, 0, 0, 0, 2, 0, 0,
//                       0, 0, 0 , 2, 0, 2, 0, 0, 1, 0,
//                       2, 0, 0 , 0, 2, 0, 0, 1, 0, 1});


  std::vector<double> big_matrix (n * n);
  srand (1);
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        {
          if (i == j)
            {
              big_matrix[i * n + j] = 1;
              continue;
            }
          if (rand  () % 4 <= 2)
            {
              big_matrix[i * n + j] = 0;
              continue;
            }
          big_matrix[i * n + j] = 2;
        }
    }
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      big_matrix[i * n + j] = big_matrix[j * n + i];

  msr.convert (n, big_matrix);

  msr_matrix save;
  save = msr;
 // msr.dump ();

  std::vector<double> x (n);

  for (int i = 0; i < n; i++)
    x[i] = (i & 1) * 10;

  msr_matrix precond;

  pthread_barrier_t barrier;
  pthread_barrier_init (&barrier, NULL, p);
  pthread_t pt;
  std::vector<double> buf (n);
  std::vector<double> v1 (n, 0), v2 (n), v3 (n), p_sized (p);
  std::vector<double> rhs (n);
  msr.mult_vector (x, rhs);
  std::vector<double> rhs_save (rhs);

  limited_deque<std::vector<double> > basis (dim);
  limited_deque<std::vector<double> > basis_derivs (dim);
  limited_deque<std::vector<double> > turns (dim);
  std::vector<double> hessenberg (dim + 2);

  std::vector<msr_thread_dqgmres_solver> handlers;
  bool flag = false;
  for (int i = 0; i < p; i++)
    {
      handlers.push_back (msr_thread_dqgmres_solver
                          (i, p, &barrier, buf, msr, precond,
                           preconditioner_type::jacobi, dim,
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

  double res = 0;
  for (int i = 0; i < n; i++)
    {
      res += pow ((rhs_comp[i] - rhs_save[i]), 2);
    }
  res = sqrt (res);
  printf ("Residual L2 : %.6le\n", res);
  return 0;
}
