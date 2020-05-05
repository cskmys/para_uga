#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <assert.h>

typedef double mat;
#define PRINT_FMT_SPECIFIER "%f "

#define NULL_PTR_CHK(x) assert(x != NULL)
#define SAME_PTR_CHK(x, y) assert(x != y)
#define ABOVE_ZERO_CHK(x) assert(x > 0.0)

#endif // UTILS_H
