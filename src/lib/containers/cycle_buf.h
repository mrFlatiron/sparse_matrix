#ifndef CYCLE_BUF_H
#define CYCLE_BUF_H

#include <vector>
#include <cstdio>

template <class T>
class cycle_buf
{
private:
  int m_capacity;
  int m_size;
  int m_pos;
  int m_vals_looped;
  int m_oldest;
  int m_newest;
  std::vector<T> m_vals;
public:
  cycle_buf (const int capacity);
  ~cycle_buf ();
  T &get_oldest ();
  T &get_newest ();
  T &get_prev ();
  T &get_next ();
  void to_start ();
  bool is_loop_done () const;
  void insert (const T &new_val);
  bool is_at_oldest () const;
  bool is_at_newest () const;
};

#endif // CYCLE_BUF_H

template <class T>
cycle_buf<T>::~cycle_buf ()
{

}

template <class T>
T &cycle_buf<T>::get_oldest ()
{
  m_pos = m_oldest;
  return m_vals[m_pos];
}

template <class T>
T &cycle_buf<T>::get_newest ()
{
  m_pos = m_newest;
  return m_vals[m_pos];
}

template <class T>
T &cycle_buf<T>::get_prev ()
{
  if (m_pos > 0)
    m_pos--;
  else
    m_pos = m_size - 1;
  return m_vals[m_pos];
}

template <class T>
T &cycle_buf<T>::get_next ()
{

  int new_pos;
  if (m_pos < m_size - 1)
    new_pos = m_pos + 1;
  else
    new_pos = 0;

  int ret_pos;

  if (!m_vals_looped)
    {
      ret_pos = m_pos;
    }
  else
    {
      ret_pos = new_pos;
      m_pos = new_pos;
    }


  m_vals_looped++;
  return m_vals[ret_pos];
}

template <class T>
void cycle_buf<T>::to_start ()
{
  m_vals_looped = 0;
  m_pos = m_oldest;
}

template <class T>
bool cycle_buf<T>::is_loop_done () const
{
  return m_vals_looped == m_size;
}

template <class T>
void cycle_buf<T>::insert (const T &new_val)
{
  if (m_size < m_capacity)
    {
      if (m_size != 0)
        m_newest++;

      m_vals[m_newest] = new_val;
      m_size++;
      return;
    }
  m_vals[m_oldest] = new_val;
  m_newest = m_oldest;
  m_oldest = (m_oldest + 1) % m_size;
}

template <class T>
bool cycle_buf<T>::is_at_oldest () const
{
  return m_pos == m_oldest;
}

template <class T>
bool cycle_buf<T>::is_at_newest () const
{
  return m_pos == m_newest;
}

template <class T>
cycle_buf<T>::cycle_buf (const int capacity) :
  m_size (0),
  m_pos (0),
  m_at_start (false),
  m_vals_looped (0),
  m_oldest (0),
  m_newest (0)
{
  if (capacity <= 0)
    {
      fprintf (stderr, "WARNING: cycle_buf::cycle_buf: size <= 0. Assuming size = 1\n", capacity);
      m_capacity = 1;
    }
  else
    m_capacity = capacity;

  m_vals.resize (m_capacity);
}
