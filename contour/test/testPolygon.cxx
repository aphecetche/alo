//
// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// @author  Laurent Aphecetche


// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// @author  Laurent Aphecetche

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <iostream>
#include "polygon.h"

using namespace o2::mch::geometry;

struct POLYGONS
{
    POLYGONS()
    {
      testPads.push_back({{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.0, 0.0}}});
      testPads.push_back({{{1.0, 3.0}, {2.0, 3.0}, {2.0, 4.0}, {1.0, 4.0}, {1.0, 3.0}}});
      testPads.push_back({{{1.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}, {1.0, 0.0}}});
      testPads.push_back({{{0.0, 1.0}, {1.0, 1.0}, {1.0, 2.0}, {0.0, 2.0}, {0.0, 1.0}}});
      testPads.push_back({{{1.0, 1.0}, {2.0, 1.0}, {2.0, 2.0}, {1.0, 2.0}, {1.0, 1.0}}});
      testPads.push_back({{{1.0, 2.0}, {2.0, 2.0}, {2.0, 3.0}, {1.0, 3.0}, {1.0, 2.0}}});
    }

    PolygonCollection<double> testPads;
    Polygon<double> polygon;
    Polygon<double> testPolygon{
      {{0.1, 0.1}, {1.1, 0.1}, {1.1, 1.1}, {2.1, 1.1}, {2.1, 3.1}, {1.1, 3.1}, {1.1, 2.1}, {0.1, 2.1}, {0.1, 0.1}}};
    Polygon<int> counterClockwisePolygon{{0, 0},
                                         {1, 0},
                                         {1, 1},
                                         {0, 1},
                                         {0, 0}};
    Polygon<int> clockwisePolygon{{0, 0},
                                  {0, 1},
                                  {1, 1},
                                  {1, 0},
                                  {0, 0}};
    Polygon<double> clockwisePolygonDouble{{0, 0},
                                           {0, 1},
                                           {1, 1},
                                           {1, 0},
                                           {0, 0}};

};

BOOST_AUTO_TEST_SUITE(o2_mch_geometry)

BOOST_FIXTURE_TEST_SUITE(polygon, POLYGONS)

BOOST_AUTO_TEST_CASE(GetYPositions)
{
  std::vector<double> xpos, ypos;

  auto p = integralPolygon(testPads, xpos, ypos);

  const std::vector<double> expected{0, 1, 2, 3, 4};
  BOOST_TEST(ypos == expected);
}

BOOST_AUTO_TEST_CASE(CreateCounterClockwiseOrientedPolygon)
{
  BOOST_CHECK(isCounterClockwiseOriented(counterClockwisePolygon));
}

BOOST_AUTO_TEST_CASE(CreateClockwiseOrientedPolygon)
{
  BOOST_CHECK(!isCounterClockwiseOriented(clockwisePolygon));
}

BOOST_AUTO_TEST_CASE(CircularTestIntegralToDoublePolygon)
{
  std::vector<double> xPositions, yPositions;
  Polygon<int> ipolygons = integralPolygon(testPolygon, xPositions, yPositions);
  Polygon<double> test = fpPolygon(ipolygons, xPositions, yPositions);
  BOOST_CHECK(test == testPolygon);
}

BOOST_AUTO_TEST_CASE(SignedArea)
{
  BOOST_CHECK_CLOSE(signedArea(testPolygon), 4.0, 0.1);
}

BOOST_AUTO_TEST_CASE(ClosingAClosedPolygonIsANop)
{
  BOOST_CHECK(testPolygon == close(testPolygon));
}

BOOST_AUTO_TEST_CASE(ClosePolygon)
{
  Polygon<int> opened{{0, 0},
                      {1, 0},
                      {1, 1},
                      {0, 1}};
  Polygon<int> expected{{0, 0},
                        {1, 0},
                        {1, 1},
                        {0, 1},
                        {0, 0}};
  auto closed = close(opened);
  BOOST_TEST(expected==closed);
}

BOOST_AUTO_TEST_CASE(ThrowIfClosingAPolygonResultInANonManhattanPolygon)
{
  Polygon<int> triangle{{0, 0},
                        {1, 0},
                        {1, 1}};

  BOOST_CHECK_THROW(close(triangle), std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

