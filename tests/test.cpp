#include "../include/mathz/Matrix.h"
#include "../include/mathz/Misc.h"
#include "../include/mathz/Vector.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

using namespace Catch;

using namespace mathz;

TEST_CASE("Radians to Degrees", "[misc]")
{
	unsigned int num_tests = 5;
	std::pair<float, float> tests[5] = {
		{ 180.f, 3.1415f },
		{ 90.f, 1.5708f },
		{ 270.f, 4.7124f },
		{ 1.f, 0.01745f },
		{ 0.f, 0.f }
	};

	for (unsigned int i = 0; i < 5; ++i)
	{
		REQUIRE(tests[i].second == Approx( radians(tests[i].first) ).margin(0.0001).epsilon(1e-10));
	}
}

TEST_CASE("Degrees to Radians", "[misc]")
{
	const unsigned int num_tests = 5;
	std::pair<float, float> tests[num_tests] = {
		{ 180.f, 3.1415f },
		{ 90.f, 1.5708f },
		{ 270.f, 4.7124f },
		{ 1.f, 0.01745f },
		{ 0.f, 0.f }
	};

	for (unsigned int i = 0; i < num_tests; ++i)
	{
		REQUIRE(tests[i].first == Approx(degrees(tests[i].second)).margin(0.01).epsilon(1e-10));
	}
}

TEST_CASE("Normalize", "[vector 4]")
{
	const unsigned int num_tests = 2;
	std::pair<Vec4, Vec4> tests[num_tests] = {
		{ Vec4(-168.027, 449.713, 88.144, -489.202), Vec4(-0.243, 0.651, 0.127, -0.708) },
		{ Vec4(-218.636, 975.642, -549.982, -814.964), Vec4(-0.156, 0.696, -0.392, -0.581) }
	};

	for (unsigned int i = 0; i < num_tests; ++i)
	{
		tests[i].first.normalize();
		REQUIRE(tests[i].first.x == Approx(tests[i].second.x).margin(0.001).epsilon(1e-10));
		REQUIRE(tests[i].first.y == Approx(tests[i].second.y).margin(0.001).epsilon(1e-10));
		REQUIRE(tests[i].first.z == Approx(tests[i].second.z).margin(0.001).epsilon(1e-10));
		REQUIRE(tests[i].first.w == Approx(tests[i].second.w).margin(0.001).epsilon(1e-10));
	}
}

TEST_CASE("Matrix Inverse", "[matrix 4x4]")
{
	Mat4 orig(
		0.1f, 0.4f, 33.f, 90.1f,
		0.15f, 3.3f, 55.f, 92.1f,
		10.1f, 1110.4f, 330.f, 900.1f,
		0.711f, 0.434f, 3763.f, 888.881f
	);

	Mat4 expected_inverse(
		-36.625839796060384856f, 39.177836385761287555f, -0.1031442088728193802f, -0.24238392214160183891f,
		0.29224725593152282103f, -0.32226779149716653452f, 0.0017522690387522774469f, 0.0019937072292893050087f,
		-0.005507383037032889493f, 0.00281126780419296924f, -0.00000647779078754478675f, 0.00027352154715688594865f,
		0.052468687208163315323f, -0.043081557819263375435f, 0.00010907081429267486132f, 0.00015998555234480827506f
	);

	Mat4 result = orig.inverse();

	for (unsigned int i = 0; i < 4; ++i)
	{
		for (unsigned int j = 0; j < 4; ++j)
		{
			REQUIRE(result[i][j] == Approx(expected_inverse[i][j]).margin(0.0001).epsilon(1e-10));
		}
	}
}

TEST_CASE("Matrix Transpose", "[matrix 4x4]")
{
	Mat4 orig(
		1.f, 3.f, 3.f, 3.f,
		2.f, 3.f, 1.f, 3.f,
		2.f, 4.f, 1.f, 3.f,
		2.f, 5.f, 1.f, 3.f
	);

	Mat4 expected_transpose(
		1.f, 2.f, 2.f, 2.f,
		3.f, 3.f, 4.f, 5.f,
		3.f, 1.f, 1.f, 1.f,
		3.f, 3.f, 3.f, 3.f
	);

	Mat4 result = orig.transpose();

	for (unsigned int i = 0; i < 4; ++i)
	{
		for (unsigned int j = 0; j < 4; ++j)
		{
			REQUIRE(result[i][j] == Approx(expected_transpose[i][j]).margin(0.0001).epsilon(1e-10));
		}
	}
}

TEST_CASE("Matrix Addition", "[matrix 4x4]")
{
	Mat4 a(
		1.f, 3.f, 3.f, 3.f,
		2.f, 3.f, 1.f, 3.f,
		2.f, 4.f, 1.f, 3.f,
		2.f, 5.f, 1.f, 3.f
	);

	Mat4 b(
		1.2f, 2.2f, 2.3f, 22.4f,
		3.2f, 3.2f, 4.4f, 5.34f,
		3.49f, 1.1f, 1.4f, 12.1f,
		77.f, 3.2f, 3.4f, 3.32f
	);

	Mat4 expected_addition(
		2.2f, 5.2f, 5.3f, 25.4f,
		5.2f, 6.2f, 5.4f, 8.34f,
		5.49f, 5.1f, 2.4f, 15.1f,
		79.f, 8.2f, 4.4f, 6.32f
	);

	Mat4 result = a + b;

	for (unsigned int i = 0; i < 4; ++i)
	{
		for (unsigned int j = 0; j < 4; ++j)
		{
			REQUIRE(result[i][j] == Approx(expected_addition[i][j]).margin(0.0001).epsilon(1e-10));
		}
	}
}

TEST_CASE("Matrix Subtraction", "[matrix 4x4]")
{
	Mat4 a(
		1.f, 3.f, 3.f, 3.f,
		2.f, 3.f, 1.f, 3.f,
		2.f, 4.f, 1.f, 3.f,
		2.f, 5.f, 1.f, 3.f
	);

	Mat4 b(
		1.2f, 2.2f, 2.3f, 22.4f,
		3.2f, 3.2f, 4.4f, 5.34f,
		3.49f, 1.1f, 1.4f, 12.1f,
		77.f, 3.2f, 3.4f, 3.32f
	);

	Mat4 expected_subtraction(
		-0.2f, 0.8f, 0.7f, -19.4f,
		-1.2f, -0.2f, -3.4f, -2.34f,
		-1.49f, 2.9f, -0.4f, -9.1f,
		-75.f, 1.8f, -2.4f, -0.32f
	);

	Mat4 result = a - b;

	for (unsigned int i = 0; i < 4; ++i)
	{
		for (unsigned int j = 0; j < 4; ++j)
		{
			REQUIRE(result[i][j] == Approx(expected_subtraction[i][j]).margin(0.0001).epsilon(1e-10));
		}
	}
}

TEST_CASE("Matrix-Matrix Multiplication", "[matrix 4x4]")
{
	Mat4 a(
		273.f, 683.f, 792.f, 865.f,
		-585.f, -281.f, 840.f, 19.f,
		-125.f, -921.f, -840.f, 284.f,
		927.f, 446.f, -966.f, -673.f
	);

	Mat4 b(
		24.f, 307.f, -34.f, 801.f,
		-894.f, 711.f, 940.f, 599.f,
		-796.f, -976.f, -758.f, 781.f,
		303.f, -451.f, 704.f, -70.f
	);

	Mat4 expected_multiplication(
		-972387.f, -593683.f, 641362.f,	1185792.f,
		-425709.f, -1207795.f, -867594.f,	17806.f,
		1575066.f, -1450.f, -24834.f, -1327724.f,
		188541.f,	1848034.f,	646158.f,	302345.f
	);

	Mat4 result = a * b;

	for (unsigned int i = 0; i < 4; ++i)
	{
		for (unsigned int j = 0; j < 4; ++j)
		{
			REQUIRE(result[i][j] == Approx(expected_multiplication[i][j]).margin(0.0001).epsilon(1e-10));
		}
	}
}

TEST_CASE("Matrix-Vector Multiplication", "[matrix 4x4]")
{
	Mat4 a(
		380.f, 881.f, 321.f, -405.f,
		-88.f, 379.f, 395.f, -722.f,
		936.f, -314.f, -846.f, 262.f,
		554.f, -970.f, 389.f, -190.f
	);

	Vec3 b(
		 265.f,
		-810.f,
		 233.f
	);

	Vec3 expected_multiplication(
		 0.725589810365f,
		-0.274276447894f,
		-0.801733441565f
	);

	Vec3 result = a * b;

	REQUIRE(result.x == expected_multiplication.x);
	REQUIRE(result.y == expected_multiplication.y);
	REQUIRE(result.z == expected_multiplication.z);
}

TEST_CASE("Matrix-Scalar Multiplication", "[matrix 4x4]")
{
	Mat4 a(
		-369.f, 434.f, 463.f, -743.f,
		-5.f, -397.f, 403.f, -822.f,
		791.f, -994.f, -760.f, -177.f,
		-233.f, -786.f, 868.f, -4
	);

	float scalar = 2.f;

	Mat4 expected_multiplication(
		-738.f, 868.f, 926.f, -1486.f,
		-10.f, -794.f, 806.f, -1644.f,
		1582.f, -1988.f, -1520.f, -354.f,
		-466.f, -1572.f, 1736, -8
	);

	Mat4 result = a * scalar;

	for (unsigned int i = 0; i < 4; ++i)
	{
		for (unsigned int j = 0; j < 4; ++j)
		{
			REQUIRE(result[i][j] == Approx(expected_multiplication[i][j]).margin(0.0001).epsilon(1e-10));
		}
	}
}