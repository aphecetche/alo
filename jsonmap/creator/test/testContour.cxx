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
#include <boost/polygon/polygon.hpp>
#include <chrono>
#include <iostream>
#include "AliMUONContour.h"
#include "AliMUONManuContourMaker.h"
#include "AliMpMotifPosition.h"
#include "AliMpSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "mapping.h"
#include "motifPosition.h"
#include "seg.h"
#include "segnumbers.h"
#include "AliMpMotifType.h"
#include "AliMpConnection.h"
#include "AliMpManuIterator.h"
#include <vector>
#include <AliMpManuIterator.h>
#include <boost/format.hpp>
#include <AliMpMotif.h>
#include "AliMUONPolygon.h"

namespace gtl = boost::polygon;
using namespace boost::polygon::operators;
using namespace gtl; // to be able to use operators on polygon vectors

typedef gtl::polygon_90_data<int> Polygon;
typedef gtl::polygon_traits<Polygon>::point_type Point;
typedef gtl::polygon_90_set_data<gtl::point_traits<Point>::coordinate_type> PolygonSet;

constexpr int NLOOP = 1;

namespace {

void dump(const Polygon& polygon)
{
  for (auto point : polygon) {
    std::cout << boost::format("%7.2f,%7.2f") % gtl::x(point) % gtl::y(point) << std::endl;
  }
}

void dump(const char* msg, const PolygonSet& polygonSet)
{
  std::cout << msg << ":" << std::endl;
  std::vector<Polygon> pols;
  polygonSet.get_polygons(pols);
  for (const auto& p: pols) {
    dump(p);
  }
}

constexpr float CONVERT{1E4};

int cm2int(double cm) {
  return static_cast<int>(std::round(cm*CONVERT));
}

}

struct CONTOURS
{

    CONTOURS()
    {
      std::vector<std::string> segnames = get_all_segmentation_names(m.ddlStore(), m.mseg());
      std::vector<std::pair<std::vector<AliMpMotifPosition*>, std::vector<AliMpMotifPosition*>>> mp = get_motifpositions(
        m.ddlStore(), m.mseg(), segnames);
      for (auto& p:mp) {
        for (auto& pos: p.first) {
          motifPositions.push_back(pos);
        }
        for (auto& pos: p.second) {
          motifPositions.push_back(pos);
        }
      }
    }

    Mapping m;
    std::vector<AliMpMotifPosition*> motifPositions;
};

PolygonSet GTLContour(AliMpSegmentation* mseg, int detElemId, int manuId)
{
  const AliMpVSegmentation* seg = mseg->GetMpSegmentationByElectronics(detElemId, manuId);
  const AliMpSlat* slat = slat_from_seg(*seg);
  assert(slat);
  AliMpMotifPosition* pos = slat->FindMotifPosition(manuId);
  AliMpVMotif* motif = pos->GetMotif();
  AliMpMotifType* motifType = motif->GetMotifType();

  PolygonSet pads;
  std::array<Point, 4> pts;
  Polygon gtlPad;

  for (Int_t i = 0; i <= 64; ++i) {
    AliMpConnection* connection = motifType->FindConnectionByGassiNum(i);

    if (connection) {
      Int_t ix = connection->GetLocalIx();
      Int_t iy = connection->GetLocalIy();

      Double_t x, y, dx, dy;

      motif->GetPadDimensionsByIndices(ix, iy, dx, dy);
      motif->PadPositionLocal(ix, iy, x, y);

      x += pos->GetPositionX() - seg->GetPositionX();
      y += pos->GetPositionY() - seg->GetPositionY();

      pts[0].x(cm2int(x - dx));
      pts[0].y(cm2int(y - dy));

      pts[1].x(cm2int(x + dx));
      pts[1].y(cm2int(y - dy));

      pts[2].x(cm2int(x + dx));
      pts[2].y(cm2int(y + dy));

      pts[3].x(cm2int(x - dx));
      pts[3].y(cm2int(y + dy));

      gtlPad.set(pts.begin(), pts.end());
      pads.insert(gtlPad);
    }
  }
  PolygonSet result;

  result |= pads;

  return result;
}

PolygonSet AliRootContour(int detElemId, int manuId)
{
  AliMUONManuContourMaker contourMaker(nullptr);
  AliMUONContour* contour = contourMaker.CreateManuContour(detElemId, manuId);
  const TObjArray* polygons = contour->Polygons();

  TIter next(polygons);
  AliMUONPolygon* poly;

  PolygonSet result;

  while ((poly = static_cast<AliMUONPolygon*>(next()))) {

    std::vector<Point> pts;
    for (int i = 0; i < poly->NumberOfVertices(); ++i) {
      pts.push_back({cm2int(poly->X(i)), cm2int(poly->Y(i))});
    }

    Polygon contour;

    contour.set(pts.begin(), pts.end());

    result.insert(contour);

  }

  return result;
}

BOOST_FIXTURE_TEST_SUITE(mch_aliroot_mapping, CONTOURS)

BOOST_AUTO_TEST_SUITE(contour)

BOOST_AUTO_TEST_CASE(NumberOfMotifPositionsIs2265)
{
  BOOST_CHECK_EQUAL(motifPositions.size(), 2265);
}

BOOST_AUTO_TEST_CASE(TimeAliRootContourCreation)
{
  std::chrono::high_resolution_clock timer;

  auto start = timer.now();

  for (int i = 0; i < NLOOP; ++i) {
    AliMpManuIterator it;
    int deId, manuId;
    while (it.Next(deId, manuId)) {
      if (deId < 500) { continue; }
      AliMUONManuContourMaker contourMaker(nullptr);
      std::unique_ptr<AliMUONContour> contour(contourMaker.CreateManuContour(deId, manuId));
    }
  }
  std::cout << "TimeAliRootContourCreation"
            << std::chrono::duration_cast<std::chrono::milliseconds>(timer.now() - start).count() << " ms"
            << std::endl;

  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(TimeGTLContourCreation)
{
  std::chrono::high_resolution_clock timer;

  auto start = timer.now();

  for (int i = 0; i < NLOOP; ++i) {
    AliMpManuIterator it;
    int deId, manuId;
    while (it.Next(deId, manuId)) {
      if (deId < 500) { continue; }
      auto p = GTLContour(m.mseg(), deId, manuId);
    }
  }
  std::cout << "TimeGTLContourCreation"
            << std::chrono::duration_cast<std::chrono::milliseconds>(timer.now() - start).count() << " ms"
            << std::endl;

  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(CompareWithBoostPolygon)
{
  AliMpManuIterator it;

  int detElemId;
  int manuId;

  while (it.Next(detElemId,manuId)) {
    if (detElemId<500) continue;
    PolygonSet custom = AliRootContour(detElemId,manuId);
    PolygonSet gtl = GTLContour(m.mseg(), detElemId,manuId);
    bool areEqual = gtl::equivalence(custom, gtl);
    if (!areEqual) {
      std::cout << "DE " << detElemId << " MANU " << manuId << std::endl;
      dump("AliRoot",custom);
      dump("GTL",gtl);
    }
    BOOST_REQUIRE(areEqual);
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
