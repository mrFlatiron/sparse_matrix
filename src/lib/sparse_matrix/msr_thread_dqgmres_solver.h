#ifndef MSR_THREAD_DQGMRES_SOLVER_H
#define MSR_THREAD_DQGMRES_SOLVER_H

#include <vector>
#include "msr_thread_handler.h"
#include "containers/limited_deque.h"

class msr_matrix;

enum class preconditioner_type
{
  identity,
  jacobi
};

enum class solver_state
{
  OK,
  TOO_SLOW,
  MAX_ITER
};

class msr_thread_dqgmres_solver : public msr_thread_handler
{
public:
  msr_matrix &m_precond; //shared precond
  preconditioner_type m_precond_type;
  int m_dim;
  int m_max_iter;
  double m_stop_criterion;
  bool &m_flag;
  std::vector<double> &m_rhs;
  limited_deque<std::vector<double>> &m_basis;
  limited_deque<std::vector<double>> &m_basis_derivs;
  limited_deque<std::vector<double>> &m_turns;
  std::vector<double> &m_hessenberg;
  std::vector<double> &m_p_sized_buf;
  std::vector<double> &m_x;
  std::vector<double> &m_v1;
  std::vector<double> &m_v2;

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
                             bool &flag,
                             std::vector<double> &rhs,
                             limited_deque<std::vector<double>> &basis_buf,
                             limited_deque<std::vector<double>> &basis_derivs_buf,
                             limited_deque<std::vector<double>> &turns_buf,
                             std::vector<double> &hessenberg_buf,
                             std::vector<double> &p_sized_buf,
                             std::vector<double> &x,
                             std::vector<double> &v1_buf,
                             std::vector<double> &v2_buf);
  ~msr_thread_dqgmres_solver ();

  std::vector<double> &p_sized_buf () const;

  msr_matrix &precond () const;

  solver_state dqgmres_solve ();
  void compute_preconditioner ();
  void apply_preconditioner ();
  int compute_hessenberg_col ();
  void apply_turn_matrices (const int cur_iter);
  void compute_turn_matrix (const int cur_iter);
  void apply_last_turn (const int cur_iter, double &g1, double &g2);
  int compute_basis_deriv (const int cur_iter);
private:
  void synchronize ();
};

#endif // MSR_THREAD_DQGMRES_SOLVER_H
