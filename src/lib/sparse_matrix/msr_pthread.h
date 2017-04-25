#ifndef MSR_PTHREAD_H
#define MSR_PTHREAD_H

#include <vector>
#include <pthread.h>
#include "threads/thread_handler.h"
#include "msr_matrix.h"
#include "containers/cycle_buf.h"

enum class preconditioner_type
{
  jacobi
};

class msr_thread_handler : public thread_handler
{
protected:
  msr_matrix &m_matrix;
  std::vector<double> &m_shared_buf;
public:
  msr_thread_handler (const int t, const int p, pthread_barrier_t *barrier,
                      std::vector<double> &shared_buf,
                      msr_matrix &matrix);
  virtual ~msr_thread_handler ();

  double aa (const int i) const;
  void aa   (const int i, const double val);
  int ja    (const int i) const;
  void ja   (const int i, const double val);

  std::vector<double> &shared_ref () const;
  msr_matrix &matrix () const;

  void mult_vector (const msr_matrix &shared_matrix, const std::vector<double> &in,
                      std::vector<double> &out /*must be resized to n*/);
};

class msr_thread_dqgmres_solver : public msr_thread_handler
{
private:
  msr_matrix &m_precond; //shared precond
  preconditioner_type m_precond_type;
  int m_dim;
  int m_max_iter;
  double m_stop_criterion;
  std::vector<double> &m_rhs;
  cycle_buf<std::vector<double>> &m_basis;
  cycle_buf<std::vector<double>> &m_basis_derivs;
  cycle_buf<std::vector<double>> &m_turns;
  std::vector<double> &m_p_sized_buf;
  std::vector<double> &m_x;
  std::vector<double> &m_v2;
  std::vector<double> &m_v3;

public:
  msr_thread_dqgmres_solver (const int t, const int p,
                             pthread_barrier_t *barrier,
                             std::vector<double> &shared_buf,
                             msr_matrix &matrix,
                             msr_matrix &precond,
                             const preconditioner_type type,
                             const int dim,
                             const int max_iter,
                             const double stop_criterion,
                             std::vector<double> &rhs,
                             cycle_buf<std::vector<double>> &basis_buf,
                             cycle_buf<std::vector<double>> &basis_derivs_buf,
                             cycle_buf<std::vector<double>> &turns_buf,
                             std::vector<double> &p_sized_buf,
                             std::vector<double> &v1_buf,
                             std::vector<double> &v2_buf,
                             std::vector<double> &v3_buf);
  ~msr_thread_dqgmres_solver ();



  msr_matrix &precond () const;

  void dqgmres_solve ();
  void compute_preconditioner ();
  void apply_preconditioner ();
private:
  void synchronize ();
};

#endif // MSR_PTHREAD_H
