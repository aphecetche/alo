#define BOOST_TEST_MODULE mch berg test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <array>
#include <string>
#include <iostream>
#include "berg.h"

BOOST_AUTO_TEST_SUITE(berg)

int refManu2ds[64] = { 62, 61, 63, 60, 59, 55, 58, 57, 56, 54, 50, 46, 42, 39, 37, 41,
                       35, 36, 33, 34, 32, 38, 43, 40, 45, 44, 47, 48, 49, 52, 51, 53,
                       7, 6, 5, 4, 2, 3, 1, 0, 9, 11, 13, 15, 17, 19, 21, 23,
                       31, 30, 29, 28, 27, 26, 25, 24, 22, 20, 18, 16, 14, 12, 10, 8 };

BOOST_AUTO_TEST_CASE(checkBergConversion)
{
  std::array<int, 64> manu2ds = getManu2Ds("bergs-run2.json", "bergs-run3.json", true);
  bool ok{ true };
  for (auto i = 0; i < manu2ds.size(); ++i) {
    bool same = (manu2ds[i] == refManu2ds[i]);
    std::cout << i << " " << manu2ds[i] << " " << refManu2ds[i] << " " << (same ? "" : "*****") << "\n";
    if (!same) {
      ok = false;
    }
  }
  BOOST_CHECK_EQUAL(true, ok);
}

BOOST_AUTO_TEST_SUITE_END()
