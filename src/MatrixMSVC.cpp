#include "VML.h"

#include <assert.h>


namespace VML
{
	Matrix::Matrix()
	{

	}

	Matrix::Matrix(const F32 pData[16])
	{
		assert((0x0000000F & (I32)pData) == 0); // Check to see if aligned to a 16

		for (I32 i = 0; i < 4; i++)
			m_cols[i] = _mm_load_ps(&pData[4 * i]);
	}

	Matrix::Matrix(const Vector & col1, const Vector & col2, const Vector & col3, const Vector & col4)
	{
		m_cols[0] = col1.m_elems;
		m_cols[1] = col2.m_elems;
		m_cols[2] = col3.m_elems;
		m_cols[3] = col4.m_elems;
	}




	std::vector<F32> Matrix::AsFloatArray()const
	{
		std::vector<F32> v;
		v.reserve(16);
		v.resize(16);

		__declspec(align(16)) F32 f[4];

		_mm_store_ps(f, m_cols[0]);
		for (int i = 0; i < 4; i++)
			v[i] = f[i];

		_mm_store_ps(f, m_cols[1]);
		for (int i = 0; i < 4; i++)
			v[4 + i] = f[i];

		_mm_store_ps(f, m_cols[2]);
		for (int i = 0; i < 4; i++)
			v[8 + i] = f[i];

		_mm_store_ps(f, m_cols[3]);
		for (int i = 0; i < 4; i++)
			v[12 + i] = f[i];

		return v;
	}




	F32 Matrix::operator()(U8 x, U8 y)const
	{
		assert(x >= 0 || x <= 3); // Check x bounds.
		assert(y >= 0 || y <= 3); // Check y bounds.

		_declspec(align(16)) float col[4];
		_mm_store_ps(col, m_cols[y]);

		return col[x];
	}





	Vector Matrix::getColumn(U8 col)const
	{
		return Vector(m_cols[col]);
	}





	Matrix Matrix::equals(const Matrix &other, F32 epsilon)const
	{
		Matrix result;

		for (I32 i = 0; i < 4; i++)
		{
			_declspec(align(16)) F32 v[4];
			m128cmp(m_cols[i], other.m_cols[i], epsilon, v);

			for (I32 j = 0; j < 4; j++)
				v[i] = (F32)(v[i] != 0);

			result.m_cols[i] = _mm_load_ps(v);
		}

		return result;
	}

	Matrix Matrix::notEquals(const Matrix &other, F32 epsilon)const
	{
		Matrix result;

		for (I32 i = 0; i < 4; i++)
		{
			_declspec(align(16)) F32 v[4];
			m128cmp(m_cols[i], other.m_cols[i], epsilon, v);

			for (I32 j = 0; j < 4; j++)
				v[i] = (F32)(v[i] == 0);

			result.m_cols[i] = _mm_load_ps(v);
		}

		return result;
	}





	Matrix Matrix::less(const Matrix &other)const
	{
		Matrix result;

		for (I32 i = 0; i < 4; i++)
		{
			_declspec(align(16)) F32 v[4];
			__m128 vResult = _mm_cmplt_ps(m_cols[i], other.m_cols[i]);

			for (I32 j = 0; j < 4; j++)
				v[i] = (F32)(v[i] != 0);

			result.m_cols[i] = _mm_load_ps(v);
		}

		return result;
	}

	Matrix Matrix::lessEq(const Matrix &other)const
	{
		Matrix result;

		for (I32 i = 0; i < 4; i++)
		{
			_declspec(align(16)) F32 v[4];
			__m128 vResult = _mm_cmple_ps(m_cols[i], other.m_cols[i]);

			for (I32 j = 0; j < 4; j++)
				v[i] = (F32)(v[i] != 0);

			result.m_cols[i] = _mm_load_ps(v);
		}

		return result;
	}

	Matrix Matrix::greater(const Matrix &other)const
	{
		Matrix result;

		for (I32 i = 0; i < 4; i++)
		{
			_declspec(align(16)) F32 v[4];
			__m128 vResult = _mm_cmpgt_ps(m_cols[i], other.m_cols[i]);

			for (I32 j = 0; j < 4; j++)
				v[i] = (F32)(v[i] != 0);

			result.m_cols[i] = _mm_load_ps(v);
		}

		return result;
	}

	Matrix Matrix::greaterEq(const Matrix &other)const
	{
		Matrix result;

		for (I32 i = 0; i < 4; i++)
		{
			_declspec(align(16)) F32 v[4];
			__m128 vResult = _mm_cmpgt_ps(m_cols[i], other.m_cols[i]);

			for (I32 j = 0; j < 4; j++)
				v[i] = (F32)(v[i] != 0);

			result.m_cols[i] = _mm_load_ps(v);
		}

		return result;
	}





	Matrix Matrix::operator+(const Matrix &rhs)const
	{
		Matrix result;

		for (I32 i = 0; i < 4; i++)
			result.m_cols[i] = _mm_add_ps(m_cols[i], rhs.m_cols[i]);

		return result;
	}

	Matrix & Matrix::operator+=(const Matrix &rhs)
	{
		for (I32 i = 0; i < 4; i++)
			m_cols[i] = _mm_add_ps(m_cols[i], rhs.m_cols[i]);

		return *this;
	}

	Matrix Matrix::operator-(const Matrix &rhs)const
	{
		Matrix result;

		for (I32 i = 0; i < 4; i++)
			result.m_cols[i] = _mm_sub_ps(m_cols[i], rhs.m_cols[i]);

		return result;
	}

	Matrix & Matrix::operator-=(const Matrix &rhs)
	{
		for (I32 i = 0; i < 4; i++)
			m_cols[i] = _mm_sub_ps(m_cols[i], rhs.m_cols[i]);

		return *this;
	}

	Matrix Matrix::operator*(const Matrix &rhs)const
	{
		Matrix result;

		for (I32 i = 0; i < 4; i++)
		{
			result.m_cols[i] = _mm_mul_ps(_MM_REPLICATE_X_PS(rhs.m_cols[i]), m_cols[0]);
			result.m_cols[i] = _MM_MADD_PS(_MM_REPLICATE_Y_PS(rhs.m_cols[i]), m_cols[1], result.m_cols[i]);
			result.m_cols[i] = _MM_MADD_PS(_MM_REPLICATE_Z_PS(rhs.m_cols[i]), m_cols[2], result.m_cols[i]);
			result.m_cols[i] = _MM_MADD_PS(_MM_REPLICATE_W_PS(rhs.m_cols[i]), m_cols[3], result.m_cols[i]);
		}

		return result;
	}

	Matrix & Matrix::operator*=(const Matrix &rhs)
	{
		*this = operator*(rhs);

		return *this;
	}

	Vector Matrix::operator*(const Vector &rhs)const
	{
		__m128 vResult;

		vResult = _mm_mul_ps(_MM_REPLICATE_X_PS(rhs.m_elems), m_cols[0]);
		vResult = _MM_MADD_PS(_MM_REPLICATE_Y_PS(rhs.m_elems), m_cols[1], vResult);
		vResult = _MM_MADD_PS(_MM_REPLICATE_Z_PS(rhs.m_elems), m_cols[2], vResult);
		vResult = _MM_MADD_PS(_MM_REPLICATE_W_PS(rhs.m_elems), m_cols[3], vResult);

		return Vector(vResult);
	}

	Matrix Matrix::operator*(F32 scalar)const
	{
		Matrix result;

		__m128 scalarBrod = _mm_set_ps1(scalar);

		for (I32 i = 0; i < 4; i++)
		{
			result.m_cols[i] = _mm_mul_ps(scalarBrod, m_cols[i]);
		}

		return result;
	}

	Matrix & Matrix::operator*=(F32 scalar)
	{
		__m128 scalarBrod = _mm_set_ps1(scalar);

		for (I32 i = 0; i < 4; i++)
		{
			m_cols[i] = _mm_mul_ps(scalarBrod, m_cols[i]);
		}

		return *this;
	}





	void Matrix::invert(Vector * pOutOptDet)
	{
		/** ftp://download.intel.com/design/PentiumIII/sml/24504301.pdf **/

		__m128 minor0, minor1, minor2, minor3;
		__m128 det;
		if (determinantVector(minor0, minor1, minor2, minor3, det))
		{
			m_cols[0] = minor0;
			m_cols[1] = minor1;
			m_cols[2] = minor2;
			m_cols[3] = minor3;

			if (pOutOptDet)
				pOutOptDet->m_elems = det;
		}
		else
			if (pOutOptDet)
				pOutOptDet->m_elems = _mm_setzero_ps();
	}

	Vector Matrix::determinant()const
	{
		__m128 m0, m1, m2, m3, det;
		if (determinantVector(m0, m1, m2, m3, det))
			return Vector(det);
		return VectorZero();
	}





	void Matrix::transpose()
	{
		_MM_TRANSPOSE4_PS(m_cols[0], m_cols[1], m_cols[2], m_cols[3]);
	}





	bool Matrix::determinantVector(__m128 & minor0, __m128 & minor1, __m128 & minor2, __m128 & minor3, __m128 & pOutDet)const
	{
		/** ftp://download.intel.com/design/PentiumIII/sml/24504301.pdf **/

		__m128 row0, row1, row2, row3;
		__m128 det, tmp1;

		row0 = m_cols[0];

		row1 = m_cols[1];
		row1 = _mm_shuffle_ps(row1, row1, _MM_SHUFFLE(1, 0, 3, 2));

		row2 = m_cols[2];

		row3 = m_cols[3];
		row3 = _mm_shuffle_ps(row3, row3, _MM_SHUFFLE(1, 0, 3, 2));

		tmp1 = _mm_mul_ps(row2, row3);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
		minor0 = _mm_mul_ps(row1, tmp1);
		minor1 = _mm_mul_ps(row0, tmp1);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
		minor0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), minor0);
		minor1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor1);
		minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);
		// -----------------------------------------------
		tmp1 = _mm_mul_ps(row1, row2);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
		minor0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor0);
		minor3 = _mm_mul_ps(row0, tmp1);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
		minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row3, tmp1));
		minor3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor3);
		minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);
		// -----------------------------------------------
		tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, 0x4E), row3);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
		row2 = _mm_shuffle_ps(row2, row2, 0x4E);
		minor0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor0);
		minor2 = _mm_mul_ps(row0, tmp1);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
		minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row2, tmp1));
		minor2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor2);
		minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);
		// -----------------------------------------------
		tmp1 = _mm_mul_ps(row0, row1);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
		minor2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor2);
		minor3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), minor3);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
		minor2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), minor2);
		minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row2, tmp1));
		// -----------------------------------------------
		tmp1 = _mm_mul_ps(row0, row3);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
		minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row2, tmp1));
		minor2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor2);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
		minor1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor1);
		minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row1, tmp1));
		// -----------------------------------------------
		tmp1 = _mm_mul_ps(row0, row2);
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
		minor1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor1);
		minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row1, tmp1));
		tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
		minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row3, tmp1));
		minor3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor3);
		// -----------------------------------------------
		det = _mm_mul_ps(row0, minor0);
		det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
		det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);

		// Check to see if det is 0.
		F32 fDet = _mm_cvtss_f32(det);
		if (fabsf(fDet) < FLOAT_EPSILON) // Equals zero... ie: Absf(det - 0) < _FLOAT_EPSILON
			return false;

		tmp1 = _mm_rcp_ss(det);
		det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
		det = _mm_shuffle_ps(det, det, 0x00);
		minor0 = _mm_mul_ps(det, minor0);
		minor1 = _mm_mul_ps(det, minor1);
		minor2 = _mm_mul_ps(det, minor2);
		minor3 = _mm_mul_ps(det, minor3);

		_MM_TRANSPOSE4_PS(minor0, minor1, minor2, minor3);

		pOutDet = det;
		return true;
	}





	Matrix MatrixIdentity()
	{
		Matrix m;

		_declspec(align(16)) F32 col0[4] = { 1.0f, 0.0f, 0.0f, 0.0f };
		_declspec(align(16)) F32 col1[4] = { 0.0f, 1.0f, 0.0f, 0.0f };
		_declspec(align(16)) F32 col2[4] = { 0.0f, 0.0f, 1.0f, 0.0f };
		_declspec(align(16)) F32 col3[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

		m.m_cols[0] = _mm_load_ps(col0);
		m.m_cols[1] = _mm_load_ps(col1);
		m.m_cols[2] = _mm_load_ps(col2);
		m.m_cols[3] = _mm_load_ps(col3);

		return m;
	}

	Matrix MatrixLookAtLH(const Vector &position, const Vector &target, const Vector &up)
	{
		Matrix mView;

		Vector zaxis = (target - position).v3Normalize();
		Vector xaxis = (up.v3Cross(zaxis)).v3Normalize();
		Vector yaxis = zaxis.v3Cross(xaxis);

		F32 m03 = -(xaxis.v4Dot(position));
		F32 m13 = -(yaxis.v4Dot(position));
		F32 m23 = -(zaxis.v4Dot(position));

		__declspec(align(16)) F32 fX[4];
		__declspec(align(16)) F32 fY[4];
		__declspec(align(16)) F32 fZ[4];

		_mm_store_ps(fX, xaxis.m_elems);
		_mm_store_ps(fY, yaxis.m_elems);
		_mm_store_ps(fZ, zaxis.m_elems);

		_declspec(align(16)) F32 col0[4] = { fX[0], fY[0], fZ[0], 0.0f };
		_declspec(align(16)) F32 col1[4] = { fX[1], fY[1], fZ[1], 0.0f };
		_declspec(align(16)) F32 col2[4] = { fX[2], fY[2], fZ[2], 0.0f };
		_declspec(align(16)) F32 col3[4] = { m03, m13, m23, 1.0f };

		mView.m_cols[0] = _mm_load_ps(col0);
		mView.m_cols[1] = _mm_load_ps(col1);
		mView.m_cols[2] = _mm_load_ps(col2);
		mView.m_cols[3] = _mm_load_ps(col3);

		return mView;
	}

	Matrix MatrixPerspectiveFOVLH(F32 fieldOfView, F32 aspectRatio, F32 near, F32 far)
	{
		Matrix mProj;

		F32 rTanFOVDiv2 = 1.0f / tanf(fieldOfView / 2.0f);
		F32 fminusn = far - near;

		_declspec(align(16)) F32 col0[4] = { (1.0f / aspectRatio) * rTanFOVDiv2, 0.0f, 0.0f, 0.0f };
		_declspec(align(16)) F32 col1[4] = { 0.0f, rTanFOVDiv2, 0.0f, 0.0f };
		_declspec(align(16)) F32 col2[4] = { 0.0f, 0.0f, far / fminusn, 1.0f };
		_declspec(align(16)) F32 col3[4] = { 0.0f, 0.0f, -(far * near) / fminusn, 0.0f };

		mProj.m_cols[0] = _mm_load_ps(col0);
		mProj.m_cols[1] = _mm_load_ps(col1);
		mProj.m_cols[2] = _mm_load_ps(col2);
		mProj.m_cols[3] = _mm_load_ps(col3);

		return mProj;
	}

	Matrix MatrixOrthoLH(F32 width, F32 height, F32 near, F32 far)
	{
		Matrix mOrtho;

		_declspec(align(16)) F32 col0[4] = { 2.0f / width, 0.0f, 0.0f, 0.0f };
		_declspec(align(16)) F32 col1[4] = { 0.0f, 2.0f / height, 0.0f, 0.0f };
		_declspec(align(16)) F32 col2[4] = { 0.0f, 0.0f, 1.0f / (far - near), near / (near - far) };
		_declspec(align(16)) F32 col3[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

		mOrtho.m_cols[0] = _mm_load_ps(col0);
		mOrtho.m_cols[1] = _mm_load_ps(col1);
		mOrtho.m_cols[2] = _mm_load_ps(col2);
		mOrtho.m_cols[3] = _mm_load_ps(col3);

		return mOrtho;
	}

	Matrix MatrixTranslation(const Vector &translation)
	{
		Matrix mTrans;

		VECTOR3F t = translation.asVector3();

		_declspec(align(16)) F32 col0[4] = { 1.0f, 0.0f, 0.0f, 0.0f };
		_declspec(align(16)) F32 col1[4] = { 0.0f, 1.0f, 0.0f, 0.0f };
		_declspec(align(16)) F32 col2[4] = { 0.0f, 0.0f, 1.0f, 0.0f };
		_declspec(align(16)) F32 col3[4] = { t.x, t.y, t.z, 1.0f };

		mTrans.m_cols[0] = _mm_load_ps(col0);
		mTrans.m_cols[1] = _mm_load_ps(col1);
		mTrans.m_cols[2] = _mm_load_ps(col2);
		mTrans.m_cols[3] = _mm_load_ps(col3);

		return mTrans;
	}

	Matrix MatrixScaling(const Vector &scaling)
	{
		Matrix mScale;

		_declspec(align(16)) F32 col0[4] = { scaling.getX(), 0.0f, 0.0f, 0.0f };
		_declspec(align(16)) F32 col1[4] = { 0.0f, scaling.getY(), 0.0f, 0.0f };
		_declspec(align(16)) F32 col2[4] = { 0.0f, 0.0f, scaling.getZ(), 0.0f };
		_declspec(align(16)) F32 col3[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

		mScale.m_cols[0] = _mm_load_ps(col0);
		mScale.m_cols[1] = _mm_load_ps(col1);
		mScale.m_cols[2] = _mm_load_ps(col2);
		mScale.m_cols[3] = _mm_load_ps(col3);

		return mScale;
	}

	Matrix MatrixRotationX(F32 angle)
	{
		Matrix mRot;

		F32 s = sinf(angle);
		F32 c = cosf(angle);

		_declspec(align(16)) F32 col0[4] = { 1.0f, 0.0f, 0.0f, 0.0f };
		_declspec(align(16)) F32 col1[4] = { 0.0f,    c,   -s, 0.0f };
		_declspec(align(16)) F32 col2[4] = { 0.0f,    s,    c, 0.0f };
		_declspec(align(16)) F32 col3[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

		mRot.m_cols[0] = _mm_load_ps(col0);
		mRot.m_cols[1] = _mm_load_ps(col1);
		mRot.m_cols[2] = _mm_load_ps(col2);
		mRot.m_cols[3] = _mm_load_ps(col3);

		return mRot;
	}

	Matrix MatrixRotationY(F32 angle)
	{
		Matrix mRot;

		F32 s = sinf(angle);
		F32 c = cosf(angle);

		_declspec(align(16)) F32 col0[4] = { c, 0.0f,    s, 0.0f };
		_declspec(align(16)) F32 col1[4] = { 0.0f, 1.0f, 0.0f, 0.0f };
		_declspec(align(16)) F32 col2[4] = { -s, 0.0f,    c, 0.0f };
		_declspec(align(16)) F32 col3[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

		mRot.m_cols[0] = _mm_load_ps(col0);
		mRot.m_cols[1] = _mm_load_ps(col1);
		mRot.m_cols[2] = _mm_load_ps(col2);
		mRot.m_cols[3] = _mm_load_ps(col3);

		return mRot;
	}

	Matrix MatrixRotationZ(F32 angle)
	{
		Matrix mRot;

		F32 s = sinf(angle);
		F32 c = cosf(angle);

		_declspec(align(16)) F32 col0[4] = { c,   -s, 0.0f, 0.0f };
		_declspec(align(16)) F32 col1[4] = { s,    c, 0.0f, 0.0f };
		_declspec(align(16)) F32 col2[4] = { 0.0f, 0.0f, 1.0f, 0.0f };
		_declspec(align(16)) F32 col3[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

		mRot.m_cols[0] = _mm_load_ps(col0);
		mRot.m_cols[1] = _mm_load_ps(col1);
		mRot.m_cols[2] = _mm_load_ps(col2);
		mRot.m_cols[3] = _mm_load_ps(col3);

		return mRot;
	}

	Matrix MatrixRotationQuaternion(const Quaternion & qRot)
	{
		Matrix m;

		F32 xy2 = 2 * qRot.v.x * qRot.v.y;
		F32 zw2 = 2 * qRot.v.z * qRot.w;
		F32 xz2 = 2 * qRot.v.x * qRot.v.z;
		F32 yw2 = 2 * qRot.v.y * qRot.w;
		F32 yz2 = 2 * qRot.v.y * qRot.v.z;
		F32 xw2 = 2 * qRot.v.x * qRot.w;

		F32 xx2 = 2 * qRot.v.x * qRot.v.x;
		F32 yy2 = 2 * qRot.v.y * qRot.v.y;
		F32 zz2 = 2 * qRot.v.z * qRot.v.z;

		VML_ALIGN_MS(16) F32 col0[4] VML_ALIGN_GCC(16) = { 1 - yy2 - zz2, xy2 - zw2, xz2 + yw2, 0.0f };
		VML_ALIGN_MS(16) F32 col1[4] VML_ALIGN_GCC(16) = { xy2 + zw2, 1 - xx2 - zz2, yz2 - xw2, 0.0f };
		VML_ALIGN_MS(16) F32 col2[4] VML_ALIGN_GCC(16) = { xz2 - yw2, yz2 + xw2, 1 - xx2 - yy2, 0.0f };
		VML_ALIGN_MS(16) F32 col3[4] VML_ALIGN_GCC(16) = { 0.0f, 0.0f, 0.0f, 1.0f };

		m.m_cols[0] = _mm_load_ps(col0);
		m.m_cols[1] = _mm_load_ps(col1);
		m.m_cols[2] = _mm_load_ps(col2);
		m.m_cols[3] = _mm_load_ps(col3);

		return m;
	}

	Matrix MatrixRotationEuler(F32 x, F32 y, F32 z)
	{
		Matrix m;

		F32 cX = cosf(x);
		F32 cY = cosf(y);
		F32 cZ = cosf(z);

		F32 sX = sinf(x);
		F32 sY = sinf(y);
		F32 sZ = sinf(z);

		VML_ALIGN_MS(16) F32 col0[4] VML_ALIGN_GCC(16) = { cY * cX, cY * sX, -sY, 0.0f };
		VML_ALIGN_MS(16) F32 col1[4] VML_ALIGN_GCC(16) = { -(cZ * sX) + (sZ * sY * cX), (cZ * cX) + (sZ * sY * sX), sZ * cY, 0.0f };
		VML_ALIGN_MS(16) F32 col2[4] VML_ALIGN_GCC(16) = { (sZ * sX) + (cZ * sY * cX), -(sZ * cY) + (cZ * sY * sX), cZ * cY, 0.0f };
		VML_ALIGN_MS(16) F32 col3[4] VML_ALIGN_GCC(16) = { 0.0f, 0.0f, 0.0f, 1.0f };

		m.m_cols[0] = _mm_load_ps(col0);
		m.m_cols[1] = _mm_load_ps(col1);
		m.m_cols[2] = _mm_load_ps(col2);
		m.m_cols[3] = _mm_load_ps(col3);

		return m;
	}
}