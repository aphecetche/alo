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
#include <boost/test/data/test_case.hpp>
#include "vertex.h"

using namespace o2::mch::contour;

BOOST_AUTO_TEST_SUITE(o2_mch_geometry)
BOOST_AUTO_TEST_SUITE(vertex)

BOOST_AUTO_TEST_CASE(Vertical)
{
  Vertex<int> v1{12, 0};
  Vertex<int> v2{12, 20};
  BOOST_TEST(isVertical(v1, v2));
  Vertex<int> v3{0, 0};
  BOOST_TEST(isVertical(v1, v3) == false);
}

BOOST_AUTO_TEST_CASE(Horizontal)
{
  Vertex<int> v1{0, 12};
  Vertex<int> v2{20, 12};
  BOOST_TEST(isHorizontal(v1, v2));
  Vertex<int> v3{0, 0};
  BOOST_TEST(isHorizontal(v1, v3) == false);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
