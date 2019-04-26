#ifndef VECTORS_H
#define VECTORS_H


#include <algorithm>

#include "bitvector.h"


namespace CSA
{


/*
  This function merges two vectors using marked positions.
  The original vectors are deleted.
  User must delete the returned vector.
*/

template<class V>
V*
mergeVectors(V* first, V* second, usint* positions, usint n, usint size, usint block_size)
{
  if((first == 0 && second == 0) || positions == 0) { return 0; }

  typename V::Iterator* first_iter = 0;
  typename V::Iterator* second_iter = 0;

  pair_type first_run;
  bool first_finished;
  if(first == 0)
  {
    first_run = pair_type(size, 0);
    first_finished = true;
  }
  else
  {
    first_iter = new typename V::Iterator(*first);
    first_run = first_iter->selectRun(0, size);
    first_run.second++;
    first_finished = false;
  }

  usint second_bit;
  if(second == 0)
  {
    second_bit = n;
  }
  else
  {
    second_iter = new typename V::Iterator(*second);
    second_bit = second_iter->select(0);
  }

  typename V::Encoder encoder(block_size);
  for(usint i = 0; i < n; i++)
  {
    while(!first_finished && first_run.first + i < positions[i])
    {
      usint bits = std::min(first_run.second, positions[i] - i - first_run.first);
      encoder.addRun(first_run.first + i, bits);
      first_run.first += bits;
      first_run.second -= bits;
      if(first_run.second == 0)
      {
        if(first_iter->hasNext())
        {
          first_run = first_iter->selectNextRun(size);
          first_run.second++;
        }
        else { first_finished = true; }
      }
    }

    if(i == second_bit) // positions[i] is one
    {
      encoder.addBit(positions[i]);
      second_bit = second_iter->selectNext();
    }
  }

  while(!first_finished)
  {
    encoder.addRun(first_run.first + n, first_run.second);
    if(first_iter->hasNext())
    {
      first_run = first_iter->selectNextRun(size);
      first_run.second++;
    }
    else { first_finished = true; }
  }

  delete first_iter; delete second_iter;
  delete first; delete second;
  encoder.flush();
  return new V(encoder, size);
}


/*
  This function encodes the ReadBuffer as the given vector type. The user must delete the returned vector.
  If there are no 1-bits in the buffer, the function returns 0.
  Does not check whether buffer is large enough.
*/

template<class V>
V*
encode(const ReadBuffer& buffer, usint block_size, usint n)
{
  usint onebits = 0;
  typename V::Encoder encoder(block_size);
  for(usint i = 0; i < n; i++)
  {
    if(buffer.isSet(i)) { encoder.addBit(i); onebits++; }
  }
  encoder.flush();

  if(onebits == 0) { return 0; }
  return new V(encoder, n);
}


} // namespace CSA


#endif // VECTORS_H
