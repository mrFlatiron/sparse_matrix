#include "thread_vector_utils.h"
#include "thread_handler.h"
#include <cmath>


void thread_utils::lin_combination_1 (thread_handler &handler,
                        std::vector<double> &shared,
                        const std::vector<double> &in,
                        std::vector<double> &out,
                        const std::vector<double> &add,
                        const double coef)
{
  const int n = in.size ();

  int begin, work;

  handler.divide_work (n, begin, work);

  for (int i = begin; i < begin + work; i++)
    {
      shared[i] = in[i] + add[i] * coef;
    }

  handler.barrier_wait ();

  for (int i = 0; i < n; i++)
    out[i] = shared[i];

  handler.barrier_wait ();
}

void thread_utils::lin_combination_1 (thread_handler &handler,
                                      std::vector<double> &shared_inout,
                                      const std::vector<double> &add,
                                      const double coef)
{
  const int n = shared_inout.size ();
  int begin, work;
  handler.divide_work (n, begin, work);

  for (int i = begin; i < begin + work; i++)
    {
      shared_inout[i] += add[i] * coef;
    }

  handler.barrier_wait ();
}

double thread_utils::l2_norm (thread_handler &handler, const std::vector<double> &vect,
                std::vector<double> &shared_buf /* p size */)
{
  int n = vect.size ();
  int begin, work;
  handler.divide_work (n, begin, work);

  shared_buf[handler.t_id ()] = 0;

  for (int i = begin; i < begin + work; i++)
    {
      shared_buf[handler.t_id ()] += vect[i] * vect[i];
    }

  handler.barrier_wait ();

  double s = 0;
  for (int i = 0; i < handler.p (); i++)
    s += shared_buf[i];

  handler.barrier_wait ();
  return sqrt (s);
}

void thread_utils::mult_vector_coef (thread_handler &handler, std::vector<double> &shared_inout,
                                     const double coef)
{
  int n = shared_inout.size ();
  int begin, work;
  handler.divide_work (n, begin, work);

  for (int i = begin; i < begin + work; i++)
    shared_inout[i] *= coef;
}

double thread_utils::dot_product (thread_handler &handler, const std::vector<double> &in1,
                                  const std::vector<double> &in2,
                                  std::vector<double> &shared_buf)
{
  int n = in1.size ();
  int begin, work;
  handler.divide_work (n, begin, work);

  double s = 0;
  for (int i = begin; i < begin + work; i++)
    s += in1[i] * in2[i];

  shared_buf[handler.t_id ()] = s;
  s = 0;

  handler.barrier_wait ();

  for (int i = 0; i < handler.p (); i++)
    s += shared_buf[i];

  handler.barrier_wait ();
  return s;
}
