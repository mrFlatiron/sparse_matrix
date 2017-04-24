#include "msr_pthread.h"
#include "msr_matrix.h"

msr_pthread::msr_pthread (msr_pthread_initializer &initer)
{
  m_matrix = initer.matrix;
  m_t = initer.t;
  m_p = initer.p;
  m_barrier = initer.barrier;
  m_shared_buf = initer.buf;
}

msr_pthread::msr_pthread (msr_matrix *matrix, const int t, const int p, pthread_barrier_t *barrier)
{
  m_matrix = matrix;
  m_t = t;
  m_p = p;
  if (m_p <= 0)
    {
      printf ("assigning m_p = 1");
      m_p = 1;
    }
  m_barrier = barrier;
}

msr_pthread::~msr_pthread()
{

}

double msr_pthread::aa (const int i) const
{
  return m_matrix->aa (i);
}

int msr_pthread::ja (const int i) const
{
  return m_matrix->ja (i);
}

int msr_pthread::t_id () const
{
  return m_t;
}

void msr_pthread::mult_vector (const std::vector<double> &in,
                               std::vector<double> &out /*must be resized to n*/)
{
  int n = m_matrix->n ();
  int work = n / m_p;
  int r = n - work * m_p;

  int i1;

  if (m_t < r)
    {
      work++;
      i1 = work * m_t;
    }
  else
    {
      i1 = work * m_t + r;
    }

  for (int i = i1; i < i1 + work; i++)
    {
      double s = 0;
      for (int ja_iter = ja (i); ja_iter < ja (i + 1); ja_iter++)
        s += aa (ja_iter) * in[ja (ja_iter)];

      s += aa (i) * in[i];
      (*m_shared_buf)[i] = s;
    }


  pthread_barrier_wait (m_barrier);

  for (int i = 0; i < n; i++)
    out[i] = (*m_shared_buf)[i];

  pthread_barrier_wait (m_barrier);
}
