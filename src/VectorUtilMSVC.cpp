#include "VML.h"

namespace VML
{
	void m128cmp(const __m128 &first, const __m128 &second, F32 epsilon, F32 pOutData[4])
	{
		__m128 vSub = _mm_sub_ps(first, second);

		__m128 vAbsResult;
		m128abs(vAbsResult);
		__m128 vEpsilon = _mm_set_ps1(epsilon);

		__m128 vResult = _mm_cmple_ps(vAbsResult, vEpsilon);
		_mm_store_ps(pOutData, vResult);
	}

	void m128abs(__m128 &vector)
	{
		static __m128 vNegZero = _mm_set_ps1(-0.0f); // 0x80000000

		vector = _mm_andnot_ps(vNegZero, vector);
	}
}