#include "msr_thread_dqgmres_solver.h"
#include "msr_matrix.h"
#include "threads/thread_vector_utils.h"
#include <cmath>

msr_thread_dqgmres_solver::msr_thread_dqgmres_solver (const int t, const int p,
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
                                                      std::vector<double> &v2_buf):
  msr_thread_handler (t, p, barrier, shared_buf, matrix),
  m_precond (precond),
  m_precond_type (type),
  m_dim (dim),
  m_max_iter (max_iter),
  m_stop_criterion (stop_criterion),
  m_flag (flag),
  m_rhs (rhs),
  m_basis (basis_buf),
  m_basis_derivs (basis_derivs_buf),
  m_turns (turns_buf),
  m_hessenberg (hessenberg_buf),
  m_p_sized_buf (p_sized_buf),
  m_x (x),
  m_v1 (v1_buf),
  m_v2 (v2_buf)
{


}

msr_thread_dqgmres_solver::~msr_thread_dqgmres_solver ()
{

}

std::vector<double> &msr_thread_dqgmres_solver::p_sized_buf () const
{
  return m_p_sized_buf;
}



msr_matrix &msr_thread_dqgmres_solver::precond () const
{
  return m_precond;
}



solver_state msr_thread_dqgmres_solver::dqgmres_solve ()
{
  const double EPS = 1e-10;

  double rhs_norm = thread_utils::l2_norm (*this, m_rhs, m_p_sized_buf);

  compute_preconditioner ();
  apply_preconditioner ();

  mult_vector_shared_out (m_precond, m_rhs, m_v2);
  thread_utils::copy_shared (*this, m_v2, m_rhs);

  double g1 = rhs_norm;


  if (g1 < EPS)
    return solver_state::OK;

  double g2;

  thread_utils::mult_vector_coef (*this, m_rhs, 1 / g1);

  if (t_id () == 0)
    m_basis.insert (m_rhs);

  barrier_wait ();

  for (int iter = 1; iter <= m_max_iter; iter++)
    {
      bool exact_solution = false;
      if (compute_hessenberg_col ())
        {
          if (t_id () == 0)
            printf ("Exact solution!\n");
          exact_solution = true;
        }
      apply_turn_matrices (iter);
      compute_turn_matrix (iter);
      apply_last_turn (iter, g1, g2);
      if (compute_basis_deriv (iter))
        {
          return solver_state::TOO_SLOW;
        }

      if (t_id () == 0)
        m_v1 = *(m_basis_derivs.get_newest ());
      barrier_wait ();
      thread_utils::lin_combination_1 (*this, m_x, m_v1, g1);
      if (exact_solution)
        return solver_state::OK;
      double residual = fabs (g2);

      if (t_id () == 0)
        printf ("\tit=%3d, residual = %le\n", iter, residual);

      if (residual < m_stop_criterion)
        return solver_state::OK;


      g1 = g2;

    }

  return solver_state::MAX_ITER;
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
      switch (m_precond_type)
        {
        case preconditioner_type::identity:
          m_precond.aa (i, 1.);
          m_precond.ja (i, n);
          break;
        case preconditioner_type::jacobi:
          m_precond.aa (i, 1 / m_matrix.aa (i));
          m_precond.ja (i, n);
          break;
        }
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
        break;
      }
    case preconditioner_type::identity:
      break;
    }
  barrier_wait ();
  return;
}

int msr_thread_dqgmres_solver::compute_hessenberg_col ()
{
  int h_iter = 1;

  mult_vector_shared_out (matrix (), *(m_basis.get_newest ()), m_v2);

  if (t_id () == 0)
    m_basis.to_preoldest ();

  while (thread_utils::limited_deque_get_next<std::vector<double> >
         (*this, m_basis, m_v1, m_flag))
    {
      double h = thread_utils::dot_product (*this, m_v1, m_v2, m_p_sized_buf);
      if (t_id () == 0)
        m_hessenberg[h_iter] = h;
      thread_utils::lin_combination_1 (*this, m_v2, m_v1, -h);
      h_iter++;
    }

  double norm = thread_utils::l2_norm (*this, m_v2, m_p_sized_buf);

  if (t_id () == 0)
    m_hessenberg[h_iter] = norm;

  if (norm < 1e-10)
    {
      if (t_id () == 0)
        m_basis.insert (m_v2);
      barrier_wait ();
      return -1;
    }

  thread_utils::mult_vector_coef (*this, m_v2, 1 / norm);

  if (t_id () == 0)
    m_basis.insert (m_v2);

  barrier_wait ();
  return 0;
}

void msr_thread_dqgmres_solver::apply_turn_matrices (const int cur_iter)
{
  if (t_id () != 0)
    {
      barrier_wait ();
      return;
    }
  int h_iter = (cur_iter <= m_dim) ? 1 : 0;
  auto turn_ptr = m_turns.get_oldest ();
  if (!turn_ptr)
    {
      barrier_wait ();
      return;
    }
  double cos = turn_ptr->at (0);
  double sin = turn_ptr->at (1);
  double x = m_hessenberg[h_iter];
  double y = m_hessenberg[h_iter + 1];

  if (h_iter == 0)
    {
      m_hessenberg[h_iter] = y * sin;
      m_hessenberg[h_iter + 1] = y * cos;
      turn_ptr = m_turns.get_next ();
      h_iter++;
    }

  for (; turn_ptr; turn_ptr = m_turns.get_next (), h_iter++)
    {
      x = m_hessenberg[h_iter];
      y = m_hessenberg[h_iter + 1];
      cos = turn_ptr->at (0);
      sin = turn_ptr->at (1);
      m_hessenberg[h_iter] = cos * x + sin * y;
      m_hessenberg[h_iter + 1] = -sin * x + cos * y;
    }

  barrier_wait ();
}

void msr_thread_dqgmres_solver::compute_turn_matrix (const int cur_iter)
{
  if (t_id () != 0)
    {
      barrier_wait ();
      return;
    }
  int h_iter = (cur_iter <= m_dim) ? cur_iter : m_dim;
  double h_mm = m_hessenberg[h_iter];
  double h_m1m = m_hessenberg[h_iter + 1];

  double s = sqrt (h_mm * h_mm + h_m1m * h_m1m);

  double cos = h_mm / s;
  double sin = h_m1m / s;

  m_turns.insert ({cos, sin});
  barrier_wait ();
}

void msr_thread_dqgmres_solver::apply_last_turn (const int cur_iter, double &g1, double &g2)
{
  if (t_id () == 0)
    {
      double cos = m_turns.get_newest ()->at (0);
      double sin = m_turns.get_newest ()->at (1);
      int h_iter = (cur_iter <= m_dim) ? cur_iter : m_dim;
      if (t_id () == 0)
        m_hessenberg[h_iter] =  cos * m_hessenberg[h_iter] + sin * m_hessenberg[h_iter + 1];

      g2 = -sin * g1;
      g1 = cos * g1;
      if (p () >= 2)
        {
          m_p_sized_buf[0] = g1;
          m_p_sized_buf[1] = g2;
        }
      barrier_wait ();
    }
  else
    {
      barrier_wait ();
      g1 = m_p_sized_buf[0];
      g2 = m_p_sized_buf[1];
    }
  barrier_wait ();
}

int msr_thread_dqgmres_solver::compute_basis_deriv (const int cur_iter)
{
  if (t_id () == 0)
    {
      m_basis.get_newest ();
      m_v2 = *(m_basis.get_prev ());
    }
    barrier_wait ();

  int h_iter = (cur_iter <= m_dim) ? 1 : 0;

  if (t_id () == 0)
    m_basis_derivs.to_preoldest ();

  while (thread_utils::limited_deque_get_next<std::vector<double> >
         (*this, m_basis_derivs, m_v1, m_flag))
    {
      thread_utils::lin_combination_1 (*this, m_v2, m_v1, -m_hessenberg[h_iter]);
      h_iter++;
    }

  printf ("h_iter = %d\n", h_iter);
  if (fabs (m_hessenberg[h_iter]) < 1e-64)
    {
      barrier_wait ();
      return -1;
    }

  thread_utils::mult_vector_coef (*this, m_v2, 1 / m_hessenberg[h_iter]);

  if (t_id () == 0)
    m_basis_derivs.insert (m_v2);

  barrier_wait ();
  return 0;
}


