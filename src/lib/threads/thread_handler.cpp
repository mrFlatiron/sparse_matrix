#include "thread_handler.h"
#include <cstdio>

thread_handler::thread_handler (const int t_id, const int p,
                          pthread_barrier_t *barrier) :
  m_t (t_id),
  m_p (p),
  m_barrier (barrier)
{
  if (m_p <= 0)
    {
      printf ("assigning m_p = 1");
      m_p = 1;
    }
}

thread_handler::~thread_handler ()
{

}

int thread_handler::t_id () const
{
  return m_t;
}

int thread_handler::p () const
{
  return m_p;
}

void thread_handler::barrier_wait ()
{
  pthread_barrier_wait (m_barrier);
}

void thread_handler::divide_work (const int n, int &begin, int &work)
{
  work = n / m_p;
  int r = n - work * m_p;

  if (m_t < r)
    {
      work++;
      begin = work * m_t;
    }
  else
    {
      begin = work * m_t + r;
    }
}
