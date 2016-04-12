#include <cmath>

int isqrt(int value)
{
  // C++11 actually has an overload for integer value, which casts value to double automatically
  int result = static_cast<int>(sqrt(static_cast<double>(value)));
  do { ++result; } while(result * result <= value);
  do { --result; } while(result * result > value);
  return result;
}