#ifndef THREAD_HANDLER_H
#define THREAD_HANDLER_H

#include <pthread.h>

class thread_handler
{
public:
  int m_t;
  int m_p;
  pthread_barrier_t *m_barrier;
public:
  thread_handler (const int t_id, const int p,
               pthread_barrier_t *barrier);
  virtual ~thread_handler ();
  int t_id  () const;
  int p     () const;
  void barrier_wait ();
  void divide_work (const int n, int &begin, int &work);
};

#endif // THREAD_HANDLER_H
