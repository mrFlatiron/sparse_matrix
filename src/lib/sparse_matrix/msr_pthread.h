#ifndef MSR_PTHREAD_H
#define MSR_PTHREAD_H

#include <vector>
#include <pthread.h>

class msr_matrix;

struct msr_pthread_initializer
{
  msr_matrix *matrix;
  int t;
  int p;
  pthread_barrier_t *barrier;
  std::vector<double> *buf;
};

enum class preconditioner_type
{
  jacobi
};

class msr_pthread
{
private:
  msr_matrix *m_matrix; //shared matrix
  int m_t; //number of this thread
  int m_p; //number of threads
  pthread_barrier_t *m_barrier;
  std::vector<double> *m_shared_buf;
public:
  explicit msr_pthread (msr_pthread_initializer &initer);
  msr_pthread (msr_matrix *matrix, const int t, const int p, pthread_barrier_t *barrier);
  ~msr_pthread ();
  double aa (const int i) const;
  int ja (const int i) const;
  int t_id () const;
  void mult_vector (const std::vector<double> &in, std::vector<double> &out /*must be resized to n*/);
  void dqgmres_solve (const preconditioner_type precond_type, double *rhs, const int max_iter,
                      const double stop_criterion);
private:
  void synchronize ();
};

#endif // MSR_PTHREAD_H
