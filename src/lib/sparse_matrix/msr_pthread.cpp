#include "msr_pthread.h"


msr_thread_handler::msr_thread_handler (const int t, const int p,
                                        pthread_barrier_t *barrier,
                                        std::vector<double> &shared_buf, msr_matrix &matrix) :
  thread_handler (t, p, barrier),
  m_matrix (matrix),
  m_shared_buf (shared_buf)
{

}

msr_thread_handler::~msr_thread_handler ()
{

}

double msr_thread_handler::aa (const int i) const
{
  return m_matrix.aa (i);
}

void msr_thread_handler::aa (const int i, const double val)
{
  m_matrix.aa (i, val);
}

int msr_thread_handler::ja (const int i) const
{
  return m_matrix.ja (i);
}

void msr_thread_handler::ja (const int i, const double val)
{
  m_matrix.ja (i, val);
}

std::vector<double> &msr_thread_handler::shared_ref () const
{
  return m_shared_buf;
}

msr_matrix &msr_thread_handler::matrix () const
{
  return m_matrix;
}

void msr_thread_handler::mult_vector (const msr_matrix &shared_matrix, const std::vector<double> &in,
                               std::vector<double> &out /*must be resized to n*/)
{
  int n = shared_matrix.n ();
  int i1, work;

  divide_work (n, i1, work);

  for (int i = i1; i < i1 + work; i++)
    {
      double s = 0;
      for (int ja_iter = shared_matrix.ja (i); ja_iter < shared_matrix.ja (i + 1); ja_iter++)
        s += shared_matrix.aa (ja_iter) * in[shared_matrix.ja (ja_iter)];

      s += shared_matrix.aa (i) * in[i];
      m_shared_buf[i] = s;
    }


  barrier_wait ();

  for (int i = 0; i < n; i++)
    out[i] = m_shared_buf[i];

  barrier_wait ();
}

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
                                                      std::vector<double> &p_sized_buf,
                                                      std::vector<double> &v1_buf,
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
  m_p_sized_buf (p_sized_buf),
  m_x (v1_buf),
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



void msr_thread_dqgmres_solver::dqgmres_solve ()
{

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
