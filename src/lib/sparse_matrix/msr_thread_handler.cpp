#include "msr_thread_handler.h"
#include "msr_matrix.h"
#include "containers/simple_vector.h"

msr_thread_handler::msr_thread_handler (const int t, const int p,
                                        pthread_barrier_t *barrier,
                                        simple_vector &shared_buf, msr_matrix &matrix) :
  thread_handler (t, p, barrier),
  m_matrix (matrix),
  m_shared_buf (shared_buf)
{

}

msr_thread_handler::~msr_thread_handler ()
{

}

int msr_thread_handler::n() const
{
  return m_matrix.n ();
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

simple_vector &msr_thread_handler::shared_ref () const
{
  return m_shared_buf;
}

msr_matrix &msr_thread_handler::matrix () const
{
  return m_matrix;
}

void msr_thread_handler::mult_vector (const msr_matrix &shared_matrix, const simple_vector &in,
                               simple_vector &out /*must be resized to n*/)
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

void msr_thread_handler::mult_vector_shared_out (const msr_matrix &shared_matrix, const simple_vector &in,
                               simple_vector &shared_out /*must be resized to n*/)
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
      shared_out[i] = s;
    }

  barrier_wait ();
}
