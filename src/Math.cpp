#include "VML.h"

#include <math.h>

namespace VML
{
	bool FEquals(F32 first, F32 second, F32 epsilon) { return (fabsf(first - second) < epsilon); }
	bool DEquals(D64 first, D64 second, D64 epsilon) { return (fabs(first - second) < epsilon); }
}