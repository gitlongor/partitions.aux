#ifndef INTEGER_ARITH
#define INTEGER_ARITH

#include <cmath>
inline int isqrt(int value)
{
  // C++11 actually has an overload for integer value, which casts value to double automatically
  int result = static_cast<int>(sqrt(static_cast<double>(value)));
  do { ++result; } while(result * result <= value);
  do { --result; } while(result * result > value);
  return result;
}

inline int idiv_ceil(int num, int den)
{// assuming num and den both being positive 
	return num / den + ((num % den >0) ? 1 : 0);
}

#endif
