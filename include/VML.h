#pragma once

#include <stdint.h>
#include <vector>

#ifdef _MSC_VER
#include <xmmintrin.h>
#endif

#ifndef VML_ALIGN_MS
#if _MSC_VER && !__INTEL_COMPILER
#define VML_ALIGN_MS(alignment) __declspec(align(alignment))
#else
#define VML_ALIGN_MS(alignment)
#endif
#endif

#ifndef VML_ALIGN_GCC
#if __GNUC__
#define VML_ALIGN_GCC(alignment) __attribute__((aligned(alignment)))
#else
#define VML_ALIGN_GCC(alignment)
#endif
#endif




#ifndef _MM_SHUFFLE_PARAM
#define _MM_SHUFFLE_PARAM(x, y, z, w) ((x) | (y << 2) | (z << 4) | (w << 6))
#endif // _MM_SHUFFLE_PARAM

#ifndef _MM_REPLICATE_X_PS
#define _MM_REPLICATE_X_PS(v) _mm_shuffle_ps((v), (v), _MM_SHUFFLE_PARAM(0, 0, 0, 0))
#endif // _MM_REPLICATE_X_PS

#ifndef _MM_REPLICATE_Y_PS
#define _MM_REPLICATE_Y_PS(v) _mm_shuffle_ps((v), (v), _MM_SHUFFLE_PARAM(1, 1, 1, 1))
#endif // _MM_REPLICATE_Y_PS

#ifndef _MM_REPLICATE_Z_PS
#define _MM_REPLICATE_Z_PS(v) _mm_shuffle_ps((v), (v), _MM_SHUFFLE_PARAM(2, 2, 2, 2))
#endif // _MM_REPLICATE_Z_PS

#ifndef _MM_REPLICATE_W_PS
#define _MM_REPLICATE_W_PS(v) _mm_shuffle_ps((v), (v), _MM_SHUFFLE_PARAM(3, 3, 3, 3))
#endif // _MM_REPLICATE_W_PS

#ifndef _MM_MADD_PS
#define _MM_MADD_PS(v1, v2, v3) _mm_add_ps(_mm_mul_ps((v1), (v2)), (v3))
#endif // _MM_MADD_PS




namespace VML
{
	typedef uint8_t byte;

	typedef int8_t I8;
	typedef uint8_t U8;

	typedef int16_t I16;
	typedef uint16_t U16;

	typedef int32_t I32;
	typedef uint32_t U32;

	typedef int64_t I64;
	typedef uint64_t U64;

	typedef float F32;
	typedef long double D64;

	const F32 PI = 3.1415927410125732421875;
	const F32 PI_DIV2 = PI / 2.0f;
	const F32 PI_DIV3 = PI / 3.0f;
	const F32 PI_DIV4 = PI / 4.0f;
	const F32 PI_DIV6 = PI / 6.0f;
	const F32 PI_MUL2 = PI * 2.0f;
	const F32 PI_3DIV2 = 3.0f * PI / 2.0f;

	const F32 FLOAT_EPSILON = 1E-5f;
	const D64 DOUBLE_EPSILON = 1E-9;
	const D64 LONGDOUBLE_EPSILON = DOUBLE_EPSILON;

	// Single precision epsilon-equality.
	// Compares first and second, to see if they equal eachother within epsilon.
	// @param first - The first number.
	// @param second - The second number.
	// @param epsilon - Epsilon.
	// @return Were the two number equal within epsilon?
	bool FEquals(F32 first, F32 second, F32 epsilon = FLOAT_EPSILON);

	// Double precision epsilon-equality.
	// Compares first and second, to see if they equal eachother within epsilon.
	// @param first - The first number.
	// @param second - The second number.
	// @param epsilon - Epsilon.
	// @return Were the two number equal within epsilon?
	bool DEquals(D64 first, D64 second, D64 epsilon = DOUBLE_EPSILON);

	/* Holds 2 components */
	template <typename T>
	struct Vector2
	{
		Vector2<T>() {}
		Vector2<T>(T x, T y) { X = x; Y = y; }
		T X;
		T Y;
	};

	/* Holds 3 components */
	template <typename T>
	struct Vector3
	{
		Vector3<T>() {}
		Vector3<T>(T x, T y, T z) { X = x; Y = y; Z = z; }
		T X;
		T Y;
		T Z;
	};

	/* Holds 4 components */
	template <typename T>
	struct Vector4
	{
		Vector4<T>() {}
		Vector4<T>(T x, T y, T z, T w) { X = x; Y = y; Z = z; W = w; }
		T X;
		T Y;
		T Z;
		T W;
	};



	




//	 __     __        _             
//	 \ \   / /__  ___| |_ ___  _ __ 
//	  \ \ / / _ \/ __| __/ _ \| '__|
//	   \ V /  __/ (__| || (_) | |   
//		\_/ \___|\___|\__\___/|_|   
//                               
	
	// 2D vector of 32 bit single floating point components.
	struct VECTOR2F
	{
		F32 x;
		F32 y;

		VECTOR2F() { }
		VECTOR2F(F32 _x, F32 _y) : x(_x), y(_y) { }
		explicit VECTOR2F(const F32 * pDataArray) : x(pDataArray[0]), y(pDataArray[1]) { }

		VECTOR2F & operator=(const VECTOR2F & other) { x = other.x; y = other.y; return *this; }
	};

	// 2D vector of 32 bit single floating point components, aligned to 16 bytes.
	VML_ALIGN_MS(16) struct VECTOR2FA : public VECTOR2F
	{
		VECTOR2FA() : VECTOR2F() { }
		VECTOR2FA(F32 _x, F32 _y) : VECTOR2F(_x, _y) { }
		explicit VECTOR2FA(const F32 * pDataArray) : VECTOR2F(pDataArray) { }

		VECTOR2FA & operator=(const VECTOR2FA & other) { x = other.x; y = other.y; return *this; }
	} VML_ALIGN_GCC(16);

	// 3D vector of 32 bit single floating point components.
	struct VECTOR3F
	{
		F32 x;
		F32 y;
		F32 z;

		VECTOR3F() { }
		VECTOR3F(F32 _x, F32 _y, F32 _z) : x(_x), y(_y), z(_z) { }
		explicit VECTOR3F(const F32 * pDataArray) : x(pDataArray[0]), y(pDataArray[1]), z(pDataArray[2]) { }

		VECTOR3F & operator=(const VECTOR3F & other) { x = other.x; y = other.y; z = other.z; return *this; }
	};

	// 3D vector of 32 bit single floating point components, aligned to 16 bytes.
	VML_ALIGN_MS(16) struct VECTOR3FA : public VECTOR3F
	{
		VECTOR3FA() : VECTOR3F() { }
		VECTOR3FA(F32 _x, F32 _y, F32 _z) : VECTOR3F(_x, _y, _z) { }
		explicit VECTOR3FA(const F32 * pDataArray) : VECTOR3F(pDataArray) { }

		VECTOR3FA & operator=(const VECTOR3FA & other) { x = other.x; y = other.y; z = other.z; return *this; }
	} VML_ALIGN_GCC(16);

	// 4D vector of 32 bit single floating point components.
	struct VECTOR4F
	{
		F32 x;
		F32 y;
		F32 z;
		F32 w;

		VECTOR4F() { }
		VECTOR4F(F32 _x, F32 _y, F32 _z, F32 _w) : x(_x), y(_y), z(_z), w(_w) { }
		explicit VECTOR4F(const F32 * pDataArray) : x(pDataArray[0]), y(pDataArray[1]), z(pDataArray[2]), w(pDataArray[3]) { }

		VECTOR4F & operator=(const VECTOR4F & other) { x = other.x; y = other.y; z = other.z; w = other.w; return *this; }
	};

	// 4D vector of 32 bit single floating point components, aligned to 16 bytes.
	VML_ALIGN_MS(16) struct VECTOR4FA : public VECTOR4F
	{
		VECTOR4FA() : VECTOR4F() { }
		VECTOR4FA(F32 _x, F32 _y, F32 _z, F32 _w) : VECTOR4F(_x, _y, _z, _w) { }
		explicit VECTOR4FA(const F32 * pDataArray) : VECTOR4F(pDataArray) { }

		VECTOR4FA & operator=(const VECTOR4FA & other) { x = other.x; y = other.y; z = other.z; w = other.w; return *this; }
	} VML_ALIGN_GCC(16);

	VML_ALIGN_MS(16) class Vector
	{
		friend class Matrix;
		friend class Quaternion;

	public:
		Vector();
		// Creates a vector with all component set equal to "all".
		explicit Vector(F32 all);
		// Creates a vector with the four components (x, y, z, w).
		Vector(F32 x, F32 y, F32 z, F32 w);
		// Creates a vector from 2D vector.
		explicit Vector(const VECTOR2F &vector);
		// Creates a vector from 2D vector.
		explicit Vector(const VECTOR2FA &vector);
		// Creates a vector from 3D vector.
		explicit Vector(const VECTOR3F &vector);
		// Creates a vector from 3D vector.
		explicit Vector(const VECTOR3FA &vector);
		// Creates a vector from 4D vector.
		explicit Vector(const VECTOR4F &vector);
		// Creates a vector from 4D vector.
		explicit Vector(const VECTOR4FA &vector);
#ifdef _MSC_VER
	private:
		explicit Vector(__m128 vector);

	public:
#endif

		// Stores the data into a VECTOR2F object.
		VECTOR2F asVector2()const;
		// Stores the data into a VECTOR2FA object.
		VECTOR2FA asVector2A()const;
		// Stores the data into a VECTOR3F object.
		VECTOR3F asVector3()const;
		// Stores the data into a VECTOR3FA object.
		VECTOR3FA asVector3A()const;
		// Stores the data into a VECTOR4F object.
		VECTOR4F asVector4()const;
		// Stores the data into a VECTOR4FA object.
		VECTOR4FA asVector4A()const;


		// Get the x component.
		F32 getX()const;
		// Get the y component.
		F32 getY()const;
		// Get the z component.
		F32 getZ()const;
		// Get the w component.
		F32 getW()const;


		// Adds the two vectors component wise.
		Vector operator+(const Vector &other)const;
		// Adds the two vectors component wise and sets the result to this.
		Vector & operator+=(const Vector &other);

		// Subtracts the two vectors component wise.
		Vector operator-(const Vector &other)const;
		// Subtracts the two vectors component wise and sets the result to this.
		Vector & operator-=(const Vector &other);

		// Multiplies the two vectors component wise.
		Vector operator*(const Vector &other)const;
		// Multiplies the two vectors component wise and sets the result to this.
		Vector operator*=(const Vector &other);

		// Multiplies each component by the scalar.
		Vector operator*(F32 scalar)const;
		// Multiplies each component by the scalar and sets the result equal to this.
		Vector operator*=(F32 scalar);

		// Divides the two vectors component wise.
		Vector operator/(const Vector &other)const;
		// Divides the two vectors component wise ands sets the result equal to this.
		Vector operator/=(const Vector &other);

		// Divides each component by the scalar.
		Vector operator/(F32 scalar)const;
		// Divides each component by the scalar and sets the result equal to this.
		Vector operator/=(F32 scalar);


		// Sets this vector to the absolute value of itself.
		// Returns a reference to this vector.
		Vector & absoluteValue();
		// Negates this vector.
		// Returns a reference to this vector.
		Vector & negate();
		// Linearly interpolates the two vectors by the amount t and returns the resulting vector.
		Vector lerp(const Vector &other, F32 t)const;



		// Returns a string representation of this vector.
		std::string toString()const;




		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 2D Vector OPERATIONS ////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// COMPARISONS ////////////////////////////////////////////////////////////////////////////

		// Checks for equality (within epsilon) of the 2D vectors. Returns a 1.0f for equals or a 0.0f for not equals in the respective component.
		VECTOR2F v2equals(const Vector &other, F32 epsilon = FLOAT_EPSILON)const;
		// Checks for equality (within epsilon) of the 2D vectors. Returns a 1.0f for not equals or a 0.0f for equals in the respective component.
		VECTOR2F v2notEquals(const Vector &other, F32 epsilon = FLOAT_EPSILON)const;

		// Compares the 2D vectors. Returns a 1.0f for a lesser value or a 0.0f otherwise in the respective component.
		VECTOR2F v2less(const Vector &other)const;
		// Compares the 2D vectors. Returns a 1.0f for a lesser or equal value or a 0.0f otherwise in the respective component.
		VECTOR2F v2lessEq(const Vector &other)const;
		// Compares the 2D vectors. Returns a 1.0f for a greater value or a 0.0f otherwise in the respective component.
		VECTOR2F v2greater(const Vector &other)const;
		// Compares the 2D vectors. Returns a 1.0f for a greater equal value or a 0.0f otherwise in the respective component.
		VECTOR2F v2greaterEq(const Vector &other)const;


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// Vector OPERATIONS //////////////////////////////////////////////////////////////////////

		// Calculates the dot product of the two 2D vectors.
		F32 v2Dot(const Vector &other)const;

		// Calculates the length of the 2D vector.
		F32 v2Length()const;
		// Calculates the squared length of the 2D vector.
		F32 v2LengthSq()const;

		// Normalizes this 2D vector.
		// Returns a reference to this vector.
		Vector & v2Normalize();

		// Calculates the angle (in radians) between the two 2D vectors. If both vectors are of unit length, use v2AngleBetweenNormals
		// for speed.
		F32 v2AngleBetween(const Vector &other)const;
		// Calculates the angle (in radians) between the two unit 2D vectors.
		F32 v2AngleBetweenNormals(const Vector &other)const;






		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 3D Vector OPERATIONS ////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// COMPARISONS ////////////////////////////////////////////////////////////////////////////

		// Checks for equality (within epsilon) of the 3D vectors. Returns a 1.0f for equals or a 0.0f for not equals in the respective component.
		VECTOR3F v3equals(const Vector &other, F32 epsilon = FLOAT_EPSILON)const;
		// Checks for equality (within epsilon) of the 3D vectors. Returns a 1.0f for not equals or a 0.0f for equals in the respective component.
		VECTOR3F v3notEquals(const Vector &other, F32 epsilon = FLOAT_EPSILON)const;

		// Compares the 3D vectors. Returns a 1.0f for a lesser value or a 0.0f otherwise in the respective component.
		VECTOR3F v3less(const Vector &other)const;
		// Compares the 3D vectors. Returns a 1.0f for a lesser or equal value or a 0.0f otherwise in the respective component.
		VECTOR3F v3lessEq(const Vector &other)const;
		// Compares the 3D vectors. Returns a 1.0f for a greater value or a 0.0f otherwise in the respective component.
		VECTOR3F v3greater(const Vector &other)const;
		// Compares the 3D vectors. Returns a 1.0f for a greater equal value or a 0.0f otherwise in the respective component.
		VECTOR3F v3greaterEq(const Vector &other)const;


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// Vector OPERATIONS //////////////////////////////////////////////////////////////////////

		// Calculates the dot product of the two 3D vectors.
		F32 v3Dot(const Vector &other)const;

		// Calculates the length of the 3D vector.
		F32 v3Length()const;
		// Calculates the squared length of the 3D vector.
		F32 v3LengthSq()const;

		// Normalizes this 3D vector.
		// Returns a reference to this vector.
		Vector & v3Normalize();

		// Calculates a vector that is perpendicular to both vectors.
		Vector v3Cross(const Vector &other)const;

		// Calculates the angle (in radians) between the two 3D vectors. If both vectors are of unit length, use v3AngleBetweenNormals
		// for speed.
		F32 v3AngleBetween(const Vector &other)const;
		// Calculates the angle (in radians) between the two unit 3D vectors.
		F32 v3AngleBetweenNormals(const Vector &other)const;






		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 4D Vector OPERATIONS ////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// COMPARISONS ////////////////////////////////////////////////////////////////////////////

		// Checks for equality (within epsilon) of the 4D vectors. Returns a 1.0f for equals or a 0.0f for not equals in the respective component.
		VECTOR4F v4equals(const Vector &other, F32 epsilon = FLOAT_EPSILON)const;
		// Checks for equality (within epsilon) of the 4D vectors. Returns a 1.0f for not equals or a 0.0f for equals in the respective component.
		VECTOR4F v4notEquals(const Vector &other, F32 epsilon = FLOAT_EPSILON)const;

		// Compares the 4D vectors. Returns a 1.0f for a lesser value or a 0.0f otherwise in the respective component.
		VECTOR4F v4less(const Vector &other)const;
		// Compares the 4D vectors. Returns a 1.0f for a lesser or equal value or a 0.0f otherwise in the respective component.
		VECTOR4F v4lessEq(const Vector &other)const;
		// Compares the 4D vectors. Returns a 1.0f for a greater value or a 0.0f otherwise in the respective component.
		VECTOR4F v4greater(const Vector &other)const;
		// Compares the 4D vectors. Returns a 1.0f for a greater equal value or a 0.0f otherwise in the respective component.
		VECTOR4F v4greaterEq(const Vector &other)const;


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// Vector OPERATIONS //////////////////////////////////////////////////////////////////////

		// Calculates the dot product of the two 4D vectors.
		F32 v4Dot(const Vector &other)const;

		// Calculates the length of the 4D vector.
		F32 v4Length()const;
		// Calculates the squared length of the 4D vector.
		F32 v4LengthSq()const;

		// Normalizes this 4D vector.
		// Returns a reference to this vector.
		Vector & v4Normalize();

		// Calculates the angle (in radians) between the two 4D vectors. If both vectors are of unit length, use v4AngleBetweenNormals
		// for speed.
		F32 v4AngleBetween(const Vector &other)const;
		// Calculates the angle (in radians) between the two unit 4D vectors.
		F32 v4AngleBetweenNormals(const Vector &other)const;


	private:
#ifdef _MSC_VER
		__m128 m_elems;
#endif

		friend Vector VectorZero();
		friend Matrix MatrixTranslation(const Vector &translation);
		friend Matrix MatrixLookAtLH(const Vector &position, const Vector &target, const Vector &up);
	} VML_ALIGN_GCC(16);

	// Returns the zero vector.
	Vector VectorZero();

	// Compares first and second for equality.
	void m128cmp(const __m128 &first, const __m128 &second, F32 epsilon, F32 pOutData[4]);

	// Computes the absolute value of the __m128 vector.
	void m128abs(__m128 &vector);







//  __  __       _        _      
// |  \/  | __ _| |_ _ __(_)_  __
// | |\/| |/ _` | __| '__| \ \/ /
// | |  | | (_| | |_| |  | |>  < 
// |_|  |_|\__,_|\__|_|  |_/_/\_\

	// A class representing a 4x4 column vector matrix.
	VML_ALIGN_MS(16) class Matrix
	{
	public:
		Matrix();
		// Creates a 4x4 matrix with the given data.
		// pData[0..3] = Column 1.
		// pData[4..7] = Column 2.
		// pData[8..11] = Column 3.
		// pData[12..15] = Column 4.
		// WARNING: pData MUST BE ALIGNED TO A 16-BYTE BOUNDARY
		explicit Matrix(const F32 pData[16]);

		/* Creates a 4x4 matrix with the 4 row vectors.
		| col1.x col2.x col3.x col4.x |
		| col1.y col2.y col3.y col4.y |
		| col1.z col2.z col3.z col4.z |
		| col1.w col2.w col3.w col4.w | */
		Matrix(const Vector & col1, const Vector & col2, const Vector & col3, const Vector & col4);

		Matrix(const Matrix & other) = default;
		Matrix & operator=(const Matrix & other) = default;

		~Matrix() = default;

		std::vector<F32> AsFloatArray()const;

		// Returns the elements at (x, y).
		F32 operator()(U8 x, U8 y)const;

		Vector getColumn(U8 col)const;

		// Checks for equality (within epsilon) in each entry of the matrix. Returns 1.0f in entry slot if equal, 0.0f otherwise.
		Matrix equals(const Matrix &other, F32 epsilon = FLOAT_EPSILON)const;
		// Checks for inequality (within epsilon) in each entry of the matrix. Returns 1.0f in entry slot if not equal, 0.0f otherwise.
		Matrix notEquals(const Matrix &other, F32 epsilon = FLOAT_EPSILON)const;

		// Compares the matrices. Returns a 1.0f for a lesser value or a 0.0f otherwise in the respective component.
		Matrix less(const Matrix &other)const;
		// Compares the matrices. Returns a 1.0f for a lesser or equal value or a 0.0f otherwise in the respective component.
		Matrix lessEq(const Matrix &other)const;
		// Compares the matrices. Returns a 1.0f for a lesser value or a 0.0f otherwise in the respective component.
		Matrix greater(const Matrix &other)const;
		// Compares the matrices. Returns a 1.0f for a lesser value or a 0.0f otherwise in the respective component.
		Matrix greaterEq(const Matrix &other)const;

		// Adds the two matrices and returns the result.
		Matrix operator+(const Matrix & rhs)const;
		// Adds the two matrices and sets the result equal to this.
		Matrix & operator+=(const Matrix & rhs);

		// Subtracts the two matrices and returns the result.
		Matrix operator-(const Matrix & rhs)const;
		// Subtracts the two matrices and sets the result equal to this.
		Matrix & operator-=(const Matrix & rhs);

		// Performs 4x4 matrix multiplication on the two matrices and returns the result.
		Matrix operator*(const Matrix &rhs)const;
		// Performs 4x4 matrix multiplication on the two matrices and sets the result equal to this.
		Matrix & operator*=(const Matrix &rhs);

		// Performs 4x4 * 4x1 matrix multplication and returns the result.
		Vector operator*(const Vector &rhs)const;

		// Multiplies each element in the matrix by the scalar value and returns the result.
		Matrix operator*(F32 scalar)const;
		// Multiplies each element in the matrix by the scalar value and sets the result equal to this.
		Matrix & operator*=(F32 scalar);

		// Attempts to invert the matrix. pOutOptDet is an optional pointer to a Vector to hold the determinant of the matrix.
		// Note: If the determinant is equal to zero, the inverse does not exist, and thus the matrix will not change.
		void invert(Vector * pOutOptDet);
		// Caculates the determinant and broadcasts the value into each component.
		Vector determinant()const;

		// Transposes this matrix.
		void transpose();

#ifdef _MSC_VER
	private:
		// Returns true if det != 0.
		bool determinantVector(__m128 & minor0, __m128 & minor1, __m128 & minor2, __m128 & minor3, __m128 & pOutDet)const;
#endif

#ifdef _MSC_VER
	private:
		__m128 m_cols[4];
#endif

		friend class Quaternion;

		friend Matrix MatrixIdentity();
		friend Matrix MatrixLookAtLH(const Vector &position, const Vector &target, const Vector &up);
		friend Matrix MatrixPerspectiveFOVLH(F32 fieldOfView, F32 aspectRatio, F32 near, F32 far);
		friend Matrix MatrixOrthoLH(F32 width, F32 height, F32 near, F32 far);
		friend Matrix MatrixTranslation(const Vector &translation);
		friend Matrix MatrixScaling(const Vector &scaling);
		friend Matrix MatrixRotationX(F32 angle);
		friend Matrix MatrixRotationY(F32 angle);
		friend Matrix MatrixRotationZ(F32 angle);
		friend Matrix MatrixRotationQuaternion(const Quaternion & qRot);
		friend Matrix MatrixRotationEuler(F32 x, F32 y, F32 z);
		friend Quaternion QuaternionFromMatrix(const Matrix &mat);
	} VML_ALIGN_GCC(16);

	// Creates a 4x4 identity matrix.
	extern Matrix MatrixIdentity();

	// Creates a left handed look at matrix.
	extern Matrix MatrixLookAtLH(const Vector &position, const Vector &target, const Vector &up);
	// Creates a left handed perspective field of view matrix.
	extern Matrix MatrixPerspectiveFOVLH(F32 fieldOfView, F32 aspectRatio, F32 near, F32 far);
	// Creates a left handed perspective field of view matrix.
	extern Matrix MatrixOrthoLH(F32 width, F32 height, F32 near, F32 far);

	// Creates a translation matrix.
	extern Matrix MatrixTranslation(const Vector &translation);
	// Creates a scaling matrix.
	extern Matrix MatrixScaling(const Vector &scaling);
	// Creates a rotation about the x-axis matrix.
	extern Matrix MatrixRotationX(F32 angle);
	// Creates a rotation about the y-axis matrix.
	extern Matrix MatrixRotationY(F32 angle);
	// Creates a rotation about the z-axis matrix.
	extern Matrix MatrixRotationZ(F32 angle);
	// Creates a rotation matrix based off the quaternion.
	extern Matrix MatrixRotationQuaternion(const Quaternion & qRot);
	// Creates a rotation matrix based off the euler angles.
	extern Matrix MatrixRotationEuler(F32 x, F32 y, F32 z);









//   ___              _                  _             
//  / _ \ _   _  __ _| |_ ___ _ __ _ __ (_) ___  _ __  
// | | | | | | |/ _` | __/ _ \ '__| '_ \| |/ _ \| '_ \ 
// | |_| | |_| | (_| | ||  __/ |  | | | | | (_) | | | |
//  \__\_\\__,_|\__,_|\__\___|_|  |_| |_|_|\___/|_| |_|
                                                     

	// Represents a 4D unit quaternion vector used for rotations.
	class Quaternion
	{
	public:
		Quaternion();
		// Creates a quaternion, where xAxis, yAxis, and zAxis make up the unit vector to 
		// rotate about, and angleInRadians is the amount to rotate.
		Quaternion(F32 xAxis, F32 yAxis, F32 zAxis, F32 angleInRadians);
		// Creates a quaternion, where axis is the unit vector to rotate about, and
		// angleInRadians is the amount to rotate.
		Quaternion(const VECTOR3F * pAxis, F32 angleInRadians);

		VECTOR3F ToEuler()const;

		// Inverts this quaternion.
		// Returns a reference to this.
		Quaternion & invert();

		// Multiplies the two quaternions and returns the result.
		Quaternion operator*(const Quaternion &other)const;
		// Multiplies the two quaternions and sets the result equal to this.
		// Returns a reference to this.
		Quaternion & operator*=(const Quaternion &other);

		// Rotates a vector about the quaternion and returns the resulting rotated vector.
		VECTOR3F operator*(const VECTOR3F &vector)const;

		// Linearly interpolates two quaternions by the amount t and returns the normalized quaternion.
		Quaternion lerp(const Quaternion &other, F32 t)const;
		// Spherically linearly interpolates tw quaternions by the amount t and returns the quaternion.
		Quaternion slerp(const Quaternion &other, F32 t)const;

		F32 getX()const;
		F32 getY()const;
		F32 getZ()const;
		F32 getW()const;

	private:
		Quaternion interpolate(const Quaternion & q1, const Quaternion & q2, F32 weightA, F32 weightB)const;

	private:
		VECTOR3F v;
		F32 w;

		friend Matrix MatrixRotationQuaternion(const Quaternion & qRot);
		friend Quaternion QuaternionFromMatrix(const Matrix &mat);
	};

	extern Quaternion QuaternionFromEuler(F32 x, F32 y, F32 z);
	extern Quaternion QuaternionFromMatrix(const Matrix &mat);
}