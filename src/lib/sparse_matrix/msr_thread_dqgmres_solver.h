#ifndef MSR_THREAD_DQGMRES_SOLVER_H
#define MSR_THREAD_DQGMRES_SOLVER_H

#include <vector>
#include "msr_thread_handler.h"
#include "containers/cycle_buf.h"

class msr_matrix;

enum class preconditioner_type
{
  jacobi
};

enum class solver_state
{
  OK,
  MAX_ITER
};

class msr_thread_dqgmres_solver : public msr_thread_handler
{
private:
  msr_matrix &m_precond; //shared precond
  preconditioner_type m_precond_type;
  int m_dim;
  int m_max_iter;
  double m_stop_criterion;
  bool &m_flag;
  std::vector<double> &m_rhs;
  cycle_buf<std::vector<double>> &m_basis;
  cycle_buf<std::vector<double>> &m_basis_derivs;
  cycle_buf<std::vector<double>> &m_turns;
  std::vector<double> &m_hessenberg;
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
                             bool &flag,
                             std::vector<double> &rhs,
                             cycle_buf<std::vector<double>> &basis_buf,
                             cycle_buf<std::vector<double>> &basis_derivs_buf,
                             cycle_buf<std::vector<double>> &turns_buf,
                             std::vector<double> &hessenberg_buf,
                             std::vector<double> &p_sized_buf,
                             std::vector<double> &x,
                             std::vector<double> &v2_buf,
                             std::vector<double> &v3_buf);
  ~msr_thread_dqgmres_solver ();



  msr_matrix &precond () const;

  solver_state dqgmres_solve ();
  void compute_preconditioner ();
  void apply_preconditioner ();
  int compute_hessenberg_col (const int cur_iter);
  void apply_turn_matrices (const int cur_iter);
private:
  void synchronize ();
};

#endif // MSR_THREAD_DQGMRES_SOLVER_H
