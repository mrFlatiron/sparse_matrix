#include "msr_thread_dqgmres_solver.h"
#include "msr_matrix.h"
#include "threads/thread_vector_utils.h"

msr_thread_dqgmres_solver::msr_thread_dqgmres_solver (const int t, const int p,
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
                                                     std::vector<double> &hessenberg_buf,
                                                      std::vector<double> &p_sized_buf,
                                                      std::vector<double> &x,
                                                      std::vector<double> &v2_buf,
                                                      std::vector<double> &v3_buf):
  msr_thread_handler (t, p, barrier, shared_buf, matrix),
  m_precond (precond),
  m_precond_type (type),
  m_dim (dim),
  m_max_iter (max_iter),
  m_stop_criterion (stop_criterion),
  m_rhs (rhs),
  m_basis (basis_buf),
  m_basis_derivs (basis_derivs_buf),
  m_turns (turns_buf),
  m_hessenberg (hessenberg_buf),
  m_p_sized_buf (p_sized_buf),
  m_x (x),
  m_v2 (v2_buf),
  m_v3 (v3_buf)
{


}

msr_thread_dqgmres_solver::~msr_thread_dqgmres_solver ()
{

}



msr_matrix &msr_thread_dqgmres_solver::precond () const
{
  return m_precond;
}



solver_state msr_thread_dqgmres_solver::dqgmres_solve (const int cur_iter)
{
  const double EPS = 1e-15;
  std::vector<double> &r = m_rhs;
  compute_preconditioner ();
  apply_preconditioner ();
  mult_vector_shared_out (m_precond, m_x, m_v2);
  thread_utils::lin_combination_1 (*this, r, m_v2, -1);
  double g1, g2;
  g1 = thread_utils::l2_norm (*this, r, m_shared_buf);

  if (g1 < EPS)
    return solver_state::OK;

  thread_utils::mult_vector_coef (*this, r, 1 / g1);

  if (t_id () == 0)
    {
      m_basis.insert (r);
    }

  barrier_wait ();

  for (int iter = 0; iter < m_max_iter; iter++)
    {
      compute_hessenberg_col (iter);
    }
}

void msr_thread_dqgmres_solver::compute_preconditioner ()
{
  int n = m_matrix.n ();
  int begin, work;
  divide_work (n, begin, work);

  if (t_id () == 0)
    {
      m_precond.set_n (n);
      m_precond.set_arr_size (n + 1);
      m_precond.ja (n, n);
    }

  barrier_wait ();

  for (int i = begin; i < begin + work; i++)
    {
      m_precond.m_aa[i] = 1 / m_matrix.aa (i);
      m_precond.m_ja[i] = n;
    }

  barrier_wait ();
}

void msr_thread_dqgmres_solver::apply_preconditioner ()
{
  switch (m_precond_type)
    {
    case preconditioner_type::jacobi:
      {
        int begin, work;
        int n = m_matrix.n ();
        divide_work (n, begin, work);
        for (int i = begin; i < begin + work; i++)
          {
            double coef = precond ().aa (i);
            for (int ja_iter = ja (i); ja_iter < ja (i + 1); ja_iter ++)
              {
                aa (ja_iter, aa (ja_iter) * coef);
              }
            aa (i, 1);
          }
        barrier_wait ();
        return;
      }
    }
  return;
}

int msr_thread_dqgmres_solver::compute_hessenberg_col (const int cur_iter)
{
  int h_iter = 1;

  mult_vector_shared_out (matrix (), m_basis.get_newest (), m_v2);
  mult_vector_shared_out (precond (), m_v2, m_v3);


  if (t_id () == 0)
    m_basis.to_start ();

  barrier_wait ();

  while (!m_basis.loop_done ())
    {
      if (t_id () == 0)
        m_v2 = m_basis.get_next ();

      double h = thread_utils::dot_product (*this, m_v2, m_v3, m_p_sized_buf);

      if (t_id () == 0)
        m_hessenberg[h_iter] = h;

      thread_utils::lin_combination_1 (*this, m_v3, m_v2, -h);

      h_iter++;
    }

  double norm = thread_utils::l2_norm (*this, m_v3, m_p_sized_buf);

  if (t_id () == 0)
    m_hessenberg[h_iter] = norm;

  barrier_wait ();

  if (norm < 1e-64)
    return -1;

  thread_utils::mult_vector_coef (*this, m_v3, 1 / norm);

  if (t_id () == 0)
    m_basis.insert (m_v3);

  barrier_wait ();
  return 0;
}
