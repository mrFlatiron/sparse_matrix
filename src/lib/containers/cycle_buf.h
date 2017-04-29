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
  bool m_preoldest;
  int m_oldest;
  int m_newest;
  std::vector<T> m_vals;
public:
  cycle_buf (const int capacity);
  ~cycle_buf ();
  void to_preoldest ();
  T *get_oldest ();
  T *get_newest ();
  T *get_prev ();
  T *get_next ();
  void insert (const T &new_val);
};



template <class T>
cycle_buf<T>::~cycle_buf ()
{

}

template<class T>
void cycle_buf<T>::to_preoldest ()
{
  m_preoldest = true;
  m_pos = m_oldest;
}

template <class T>
T *cycle_buf<T>::get_oldest ()
{
  m_preoldest = false;
  if (m_size == 0)
    return nullptr;

  m_pos = m_oldest;
  return m_vals.data () + m_pos;
}

template <class T>
T *cycle_buf<T>::get_newest ()
{
  m_preoldest = false;
  if (m_size == 0)
    return nullptr;

  m_pos = m_newest;
  return m_vals.data () + m_pos;
}

template <class T>
T *cycle_buf<T>::get_prev ()
{
  m_preoldest = false;
  if  (m_pos == m_oldest || m_size == 0)
    return nullptr;

  if (m_pos > 0)
    m_pos--;
  else
    m_pos = m_size - 1;

  return m_vals.data () + m_pos;
}

template <class T>
T *cycle_buf<T>::get_next ()
{
  if ((m_pos == m_newest && !m_preoldest) || m_size == 0)
    return nullptr;

  if (m_preoldest)
    {
      m_pos = m_oldest;
      m_preoldest = false;
      return m_vals.data () + m_pos;
    }

  if (m_pos == m_size - 1)
    m_pos = 0;
  else
    m_pos++;

  return m_vals.data () + m_pos;
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
cycle_buf<T>::cycle_buf (const int capacity) :
  m_size (0),
  m_pos (0),
  m_preoldest (false),
  m_oldest (0),
  m_newest (0)
{
  if (capacity <= 0)
    {
      fprintf (stderr, "WARNING: cycle_buf::cycle_buf: size <= 0. Assuming size = 1\n");
      m_capacity = 1;
    }
  else
    m_capacity = capacity;

  m_vals.resize (m_capacity);
}
#endif // CYCLE_BUF_H
