#include "msr_matrix.h"
#include <cmath>

msr_matrix::msr_matrix ()
{

}

msr_matrix::~msr_matrix()
{

}

void msr_matrix::dump (FILE *fout)
{
  fprintf (fout, "msr_matrix::dump ()\n");
  fprintf (fout, "===================\n");
  fprintf (fout, "Matrix size: %d\n", m_n);
  fprintf (fout, "Arrays' size: %d\n", m_arr_size);

  fprintf (fout, "AA: ");
  for (int i = 0; i < m_arr_size; i++)
    fprintf (fout, "%.3lf ", m_aa[i]);

  fprintf (fout, "\nJA: ");

  for (int i = 0; i< m_arr_size; i++)
    fprintf (fout, "%d ", m_ja[i]);

  fprintf (fout, "\nConventional view:\n");

  for (int i = 0; i < m_n; i++)
    {
      int row_begin, row_end;
      get_ja_row_bounds (i, row_begin, row_end);
      print_row (fout, i , row_begin, row_end);
    }
  fprintf (fout, "===================\n");
}

void msr_matrix::convert (const int n, std::vector<double> matrix)
{
  m_n = n;
  int ja_iter = m_n + 1;
  int non_zero = 0;

  for (int i = 0; i < n * n; i++)
    {
      if (fabs (matrix[i]) > 1e-16)
        non_zero++;
    }

  m_arr_size = non_zero + 1;
  m_aa.resize (m_arr_size);
  m_ja.resize (m_arr_size);

  for (int i = 0; i < n; i++)
    {
      if (fabs (matrix[i * n + i] < 1e-16))
        {
          printf ("ERROR: zero diagonal element at i = %d\n", i);
          return;
        }
      m_aa[i]  = matrix[i * n + i];
    }
  m_ja[m_n] = m_arr_size;

  for (int i = 0; i < n; i++)
    {
      bool empty = true;
      for (int j = 0; j < i; j++)
        {
          if (fabs (matrix[i * n + j]) > 1e-16)
            {
              m_aa[ja_iter] = matrix[i *n + j];
              m_ja[ja_iter] =  j;
              if (empty)
                {
                  m_ja[i] = ja_iter;
                  empty = false;
                }
              ja_iter++;
            }
        }

      for (int j = i + 1; j < n; j++)
        {
          if (fabs (matrix[i * n + j]) > 1e-16)
            {
              m_aa[ja_iter] = matrix[i * n + j];
              m_ja[ja_iter] = j;
              if (empty)
                {
                  m_ja[i] = ja_iter;
                  empty = false;
                }
              ja_iter++;
            }
        }
      if (empty)
        m_ja[i] = m_arr_size;
      else
        {
          for (int k = i - 1; k >= 0; k--)
            {
              if (m_ja[k] != m_arr_size)
                break;
              m_ja[k] = m_ja[i];
            }
        }
    }
}

int msr_matrix::n () const
{
  return m_n;
}

int msr_matrix::arr_size () const
{
  return m_arr_size;
}

double msr_matrix::aa (const int i) const
{
  return m_aa[i];
}

int msr_matrix::ja (const int i) const
{
  return m_ja[i];
}


void msr_matrix::print_row (FILE *fout, const int i, const int row_begin, const int row_end)
{
  int ja_iter = row_begin;
  for (int j = 0; j < i; j++)
    {
      if (j == m_ja[ja_iter] && ja_iter < row_end)
        {
          fprintf (fout, "%.3lf ", m_aa[ja_iter]);
          ja_iter++;
          continue;
        }
      fprintf (fout, "0.000 ");
    }
  fprintf (fout, "%.3lf ", m_aa[i]);
  for (int j = i + 1; j < m_n; j++)
    {
      if (j == m_ja[ja_iter] && ja_iter < row_end)
        {
          fprintf (fout, "%.3lf ", m_aa[ja_iter]);
          ja_iter++;
          continue;
        }
      fprintf (fout, "0.000 ");
    }
  fprintf (fout, "\n");
}

void msr_matrix::get_ja_row_bounds (const int i, int &begin, int &end)
{
  begin = m_ja[i];
  end = m_ja[i + 1];
}

