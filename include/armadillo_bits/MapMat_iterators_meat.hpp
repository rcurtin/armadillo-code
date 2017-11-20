// Copyright 2008-2017 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2017 National ICT Australia (NICTA)
// Copyright 2017 Ryan Curtin (ryan@ratml.org)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup MapMat
//! @{



template<typename eT>
inline
MapMat<eT>::const_iterator::const_iterator() :
  m(NULL),
  internal_pos(0)
  {
  }



template<typename eT>
inline
MapMat<eT>::const_iterator::const_iterator(const MapMat<eT>& in, const bool end) :
  m(&in),
  internal_it(end ? m->map_ptr->end() : m->map_ptr->begin()),
  internal_pos(end ? m->map_ptr->size() : 0)
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
inline
MapMat<eT>::const_iterator::const_iterator(const MapMat<eT>& in,
                                           const uword row,
                                           const uword col) :
  m(&in),
  internal_it(m->map_ptr->begin()),
  internal_pos(0)
  {
  arma_extra_debug_sigprint_this(this);

  // Advance the position until it is in the right place.
  while (this->col() < col || (this->col() == col && this->row() < row))
    {
    ++internal_it;
    ++internal_pos;
    }
  }



template<typename eT>
inline
MapMat<eT>::const_iterator::const_iterator(const iterator& in) :
  const_iterator(*in.m, in.row(), in.col())
  {
  }



template<typename eT>
inline
typename MapMat<eT>::const_iterator&
MapMat<eT>::const_iterator::operator++()
  {
  ++internal_it;
  ++internal_pos;
  return *this;
  }



template<typename eT>
inline
typename MapMat<eT>::const_iterator&
MapMat<eT>::const_iterator::operator--()
  {
  --internal_it;
  --internal_pos;
  return *this;
  }



template<typename eT>
inline
typename MapMat<eT>::const_iterator
MapMat<eT>::const_iterator::operator++(int)
  {
  const_iterator i = *this;
  ++internal_it;
  ++internal_pos;
  return i;
  }



template<typename eT>
inline
typename MapMat<eT>::const_iterator
MapMat<eT>::const_iterator::operator--(int)
  {
  const_iterator i = *this;
  --internal_it;
  --internal_pos;
  return i;
  }



template<typename eT>
inline
eT
MapMat<eT>::const_iterator::operator*() const
  {
  return (*internal_it).second;
  }



template<typename eT>
inline
uword
MapMat<eT>::const_iterator::row() const
  {
  if (internal_it == m->map_ptr->end())
    return m->n_rows;
  else
    return (*internal_it).first % m->n_rows;
  }



template<typename eT>
inline
uword
MapMat<eT>::const_iterator::col() const
  {
  if (internal_it == m->map_ptr->end())
    return m->n_cols;
  else
    return (*internal_it).first / m->n_rows;
  }



template<typename eT>
inline
uword
MapMat<eT>::const_iterator::pos() const
  {
  return internal_pos;
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::const_iterator::operator==(const const_iterator& other) const
  {
  return (m == other.m) ? (internal_it == other.internal_it) :
      (row() == other.row() && col() == other.col());
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::const_iterator::operator!=(const const_iterator& other) const
  {
  return !(*this == other);
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::const_iterator::operator==(const typename SpSubview<eT>::const_iterator& other) const
  {
  return (col() == other.col() &&
          row() == other.row());
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::const_iterator::operator!=(const typename SpSubview<eT>::const_iterator& other) const
  {
  return !(*this == other);
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::const_iterator::operator==(const typename SpSubview<eT>::iterator& other) const
  {
  return (col() == other.col() &&
          row() == other.row());
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::const_iterator::operator!=(const typename SpSubview<eT>::iterator& other) const
  {
  return !(*this == other);
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::const_iterator::operator==(const iterator& other) const
  {
  return (m == other.m) ? (internal_it == other.internal_it) :
      (row() == other.row() && col() == other.col());
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::const_iterator::operator!=(const iterator& other) const
  {
  return !(*this == other);
  }



template<typename eT>
inline
MapMat<eT>::iterator::iterator() :
  m(NULL),
  internal_pos(0)
  {
  }



template<typename eT>
inline
MapMat<eT>::iterator::iterator(MapMat<eT>& in, const bool end) :
  m(&in),
  internal_it(end ? m->map_ptr->end() : m->map_ptr->begin()),
  internal_pos(end ? m->map_ptr->size() : 0)
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
inline
MapMat<eT>::iterator::iterator(MapMat<eT>& in,
                               const uword row,
                               const uword col) :
  m(&in),
  internal_it(m->map_ptr->begin()),
  internal_pos(0)
  {
  arma_extra_debug_sigprint_this(this);

  // Advance the position until it is in the right place.
  while (this->col() < col || (this->col() == col && this->row() < row))
    {
    ++internal_it;
    ++internal_pos;
    }
  }



template<typename eT>
inline
typename MapMat<eT>::iterator&
MapMat<eT>::iterator::operator++()
  {
  ++internal_it;
  ++internal_pos;
  return *this;
  }



template<typename eT>
inline
typename MapMat<eT>::iterator&
MapMat<eT>::iterator::operator--()
  {
  --internal_it;
  --internal_pos;
  return *this;
  }



template<typename eT>
inline
typename MapMat<eT>::iterator
MapMat<eT>::iterator::operator++(int)
  {
  iterator i = *this;
  ++internal_it;
  ++internal_pos;
  return i;
  }



template<typename eT>
inline
typename MapMat<eT>::iterator
MapMat<eT>::iterator::operator--(int)
  {
  iterator i = *this;
  --internal_it;
  --internal_pos;
  return i;
  }



template<typename eT>
inline
MapMat_val<eT>
MapMat<eT>::iterator::operator*()
  {
  return MapMat_val<eT>(*m, (*internal_it).first);
  }



template<typename eT>
inline
uword
MapMat<eT>::iterator::row() const
  {
  if (internal_it == m->map_ptr->end())
    return m->n_rows;
  else
    return (*internal_it).first % m->n_rows;
  }



template<typename eT>
inline
uword
MapMat<eT>::iterator::col() const
  {
  if (internal_it == m->map_ptr->end())
    return m->n_cols;
  else
    return (m->n_rows == 0) ? 0 : (*internal_it).first / m->n_rows;
  }



template<typename eT>
inline
uword
MapMat<eT>::iterator::pos() const
  {
  return internal_pos;
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::iterator::operator==(const const_iterator& other) const
  {
  return (m == other.m) ? (internal_it == other.internal_it) :
      (row() == other.row() && col() == other.col());
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::iterator::operator!=(const const_iterator& other) const
  {
  return !(*this == other);
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::iterator::operator==(const typename SpSubview<eT>::const_iterator& other) const
  {
  return (col() == other.col() &&
          row() == other.row());
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::iterator::operator!=(const typename SpSubview<eT>::const_iterator& other) const
  {
  return !(*this == other);
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::iterator::operator==(const typename SpSubview<eT>::iterator& other) const
  {
  return (col() == other.col() &&
          row() == other.row());
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::iterator::operator!=(const typename SpSubview<eT>::iterator& other) const
  {
  return !(*this == other);
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::iterator::operator==(const iterator& other) const
  {
  return (m == other.m) ? (internal_it == other.internal_it) :
      (row() == other.row() && col() == other.col());
  }



template<typename eT>
arma_inline
bool
MapMat<eT>::iterator::operator!=(const iterator& other) const
  {
  return !(*this == other);
  }
//! @}
