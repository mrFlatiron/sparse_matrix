#ifndef MSR_MATRIX_H
#define MSR_MATRIX_H

#include "containers/simple_vector.h"
#include <cstdio>
#include <vector>

class msr_thread_dqgmres_solver;

class msr_matrix
{
private:
  int m_n;
  int m_arr_size;
  simple_vector m_aa;
  std::vector<int> m_ja;
public:
  friend class msr_thread_dqgmres_solver;
  msr_matrix ();
  ~msr_matrix ();
  void dump (FILE *fout = stdout);
  void convert (const int n, simple_vector matrix);
  int n () const;
  void set_n (const int n);
  int arr_size () const;
  void set_arr_size (const int size);
  double aa (const int i) const;
  void aa (const int i, const double val);
  int ja (const int i) const;
  void ja (const int i, const double val);
  void set_diagonal (const simple_vector &diag_vals);
  void mult_vector (const simple_vector &in, simple_vector &out);
private:
  void print_row (FILE *fout, const int i, const int row_begin, const int row_end);
  void get_ja_row_bounds (const int i, int &begin, int &end);
};

#endif // MSR_MATRIX_H
