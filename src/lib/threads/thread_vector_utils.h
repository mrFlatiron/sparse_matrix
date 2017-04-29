#ifndef THREAD_VECTOR_UTILS_H
#define THREAD_VECTOR_UTILS_H

#include <vector>

class thread_handler;

namespace thread_utils
{
//out = in + coef * add;
void lin_combination_1 (thread_handler &handler, std::vector<double> &shared_buf,
                        const std::vector<double> &in, std::vector<double> &out,
                        const std::vector<double> &add, const double coef);

void lin_combination_1 (thread_handler &handler,
                                      std::vector<double> &shared_inout,
                                      const std::vector<double> &add,
                                      const double coef);

double l2_norm (thread_handler &handler, const std::vector<double> &vect,
                std::vector<double> &shared_buf /* p size */);

double dot_product (thread_handler &handler, const std::vector<double> &in1,
                    const std::vector<double> &in2,
                    std::vector<double> &shared_buf /* p size */);

void mult_vector_coef (thread_handler &handler, std::vector<double> &shared_inout,
                       const double coef);

void copy_shared (thread_handler &handler, const std::vector<double> &shared_source,
                  std::vector<double> &shared_out);
}

#endif // THREAD_VECTOR_UTILS_H
