#include "VML.h"

namespace VML
{
	inline F32 SIGN(float x) { return (x >= 0.0f) ? +1.0f : -1.0f; }
	inline F32 NORM(float a, float b, float c, float d) { return sqrt(a * a + b * b + c * c + d * d); }

	Quaternion QuaternionFromMatrix(const Matrix &rotMatrix)
	{
		Quaternion quat;
		__declspec(align(16)) F32 R[4][4];
		for (U8 i = 0; i < 4; i++)
			_mm_store_ps(R[i], rotMatrix.m_cols[i]);

		F32 q[4];

		q[0] = (R[0][0] + R[1][1] + R[2][2] + 1.0f) / 4.0f;
		q[1] = (R[0][0] - R[1][1] - R[2][2] + 1.0f) / 4.0f;
		q[2] = (-R[0][0] + R[1][1] - R[2][2] + 1.0f) / 4.0f;
		q[3] = (-R[0][0] - R[1][1] + R[2][2] + 1.0f) / 4.0f;
		if (q[0] < 0.0f) q[0] = 0.0f;
		if (q[1] < 0.0f) q[1] = 0.0f;
		if (q[2] < 0.0f) q[2] = 0.0f;
		if (q[3] < 0.0f) q[3] = 0.0f;
		q[0] = sqrt(q[0]);
		q[1] = sqrt(q[1]);
		q[2] = sqrt(q[2]);
		q[3] = sqrt(q[3]);
		if (q[0] >= q[1] && q[0] >= q[2] && q[0] >= q[3]) {
			q[0] *= +1.0f;
			q[1] *= SIGN(R[2][1] - R[1][2]);
			q[2] *= SIGN(R[0][2] - R[2][0]);
			q[3] *= SIGN(R[1][0] - R[0][1]);
		}
		else if (q[1] >= q[0] && q[1] >= q[2] && q[1] >= q[3]) {
			q[0] *= SIGN(R[2][1] - R[1][2]);
			q[1] *= +1.0f;
			q[2] *= SIGN(R[1][0] + R[0][1]);
			q[3] *= SIGN(R[0][2] + R[2][0]);
		}
		else if (q[2] >= q[0] && q[2] >= q[1] && q[2] >= q[3]) {
			q[0] *= SIGN(R[0][2] - R[2][0]);
			q[1] *= SIGN(R[1][0] + R[0][1]);
			q[2] *= +1.0f;
			q[3] *= SIGN(R[2][1] + R[1][2]);
		}
		else if (q[3] >= q[0] && q[3] >= q[1] && q[3] >= q[2]) {
			q[0] *= SIGN(R[1][0] - R[0][1]);
			q[1] *= SIGN(R[2][0] + R[0][2]);
			q[2] *= SIGN(R[2][1] + R[1][2]);
			q[3] *= +1.0f;
		}

		F32 r = NORM(q[0], q[1], q[2], q[3]);
		q[0] /= r;
		q[1] /= r;
		q[2] /= r;
		q[3] /= r;

		quat.w = q[0];
		quat.v.x = q[1];
		quat.v.y = q[2];
		quat.v.z = q[3];

		return quat;
	}


	Quaternion Quaternion::interpolate(const Quaternion & q1, const Quaternion & q2, F32 weightA, F32 weightB)const
	{
		Quaternion q;

		Vector v1(q1.v.x, q1.v.y, q1.v.z, q1.w);
		Vector v2(q2.v.x, q2.v.y, q2.v.z, q2.w);

		v1 *= weightA;
		v2 *= weightB;

		v1 += v2;

		v1.v4Normalize();

		__declspec(align(16)) F32 vector[4];
		_mm_store_ps(vector, v1.m_elems);

		q.v.x = vector[0];
		q.v.y = vector[1];
		q.v.z = vector[2];
		q.w = vector[3];

		return q;
	}
}