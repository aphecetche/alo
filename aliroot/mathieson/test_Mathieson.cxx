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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "AliMUONMathieson.h"

BOOST_AUTO_TEST_CASE(TestOneMathiesonIntegration)
{
  AliMUONMathieson m;

  m.SetPitch(0.21);
  m.SetSqrtKx3AndDeriveKx2Kx4(0.700);
  m.SetSqrtKy3AndDeriveKy2Ky4(0.755);

  double expected = -42.42;
  double x1 = 0;
  double y1 = 0;
  double x2 = 1;
  double y2 = 2;
  double c = m.IntXY(x1,y1,x2,y2);

  BOOST_TEST(c==expected);
}


