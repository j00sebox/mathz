#include "../Matrix.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

using namespace Catch;

TEST_CASE("Matrix Inverse", "[matrix 4x4]")
{
	mathz::Mat4 orig(
		0.1f, 0.4f, 33.f, 90.1f,
		0.15f, 3.3f, 55.f, 92.1f,
		10.1f, 1110.4f, 330.f, 900.1f,
		0.711f, 0.434f, 3763.f, 888.881f
	);

	mathz::Mat4 expected_inverse(
		-36.625839796060384856f, 39.177836385761287555f, -0.1031442088728193802f, -0.24238392214160183891f,
		0.29224725593152282103f, -0.32226779149716653452f, 0.0017522690387522774469f, 0.0019937072292893050087f,
		-0.005507383037032889493f, 0.00281126780419296924f, -0.00000647779078754478675f, 0.00027352154715688594865f,
		0.052468687208163315323f, -0.043081557819263375435f, 0.00010907081429267486132f, 0.00015998555234480827506f
	);

	mathz::Mat4 result = orig.inverse();

	for (unsigned int i = 0; i < 4; ++i)
	{
		for (unsigned int j = 0; j < 4; ++j)
		{
			REQUIRE(result(i, j) == Approx(expected_inverse(i, j)).margin(0.0001).epsilon(1e-10));
		}
	}
}