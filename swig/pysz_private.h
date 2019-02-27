#pragma once
#include "sz.h"

template <class T>
class SZTypeToTypeID
{};

#define MAKE_TYPE_TO_ID(type, id)                                              \
  template <>                                                                  \
  class SZTypeToTypeID<type>                                                   \
  {                                                                            \
  public:                                                                      \
    static const int value = id;                                               \
  };

MAKE_TYPE_TO_ID(float, SZ_FLOAT);
MAKE_TYPE_TO_ID(double, SZ_DOUBLE);

//TODO add the other types as they are supported by the dispatch method
