#ifndef MSR_MATRIX_H
#define MSR_MATRIX_H

#include <vector>
#include <cstdio>

class msr_pthread;

class msr_matrix
{
private:
  int m_n;
  int m_arr_size;
  std::vector<double> m_aa;
  std::vector<int> m_ja;
public:
//  friend class msr_pthread;
  msr_matrix ();
  ~msr_matrix ();
  void dump (FILE *fout = stdout);
  void convert (const int n, std::vector<double> matrix);
  int n () const;
  int arr_size () const;
  double aa (const int i) const;
  int ja (const int i) const;
private:
  void print_row (FILE *fout, const int i, const int row_begin, const int row_end);
  void get_ja_row_bounds (const int i, int &begin, int &end);
};

#endif // MSR_MATRIX_H
