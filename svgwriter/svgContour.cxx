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

#include "svgContour.h"
#include "contourCreator.h"

namespace o2 {
namespace mch {
namespace svg {

void writeContour(std::ofstream &out, int scale, const o2::mch::contour::Contour<double> &contour,
                  o2::mch::contour::BBox<double> &box, const char *style)
{
  for (auto i = 0; i < contour.size(); ++i) {
    out << "<polygon points=\"";
    auto vertices = getVertices(contour[i]);
    for (auto j = 0; j < vertices.size(); ++j) {
      auto v = vertices[j];
      v.x -= box.xmin();
      v.y -= box.ymin();
      out << scale * v.x << "," << scale * v.y << ' ';
    }
    out << "\" style=\"" << style << "\"/>\n";

//    if (!highlight) {
//      out << "\" style=\"fill:none;stroke:black;stroke-width:0.5px\"/>\n";
//    } else {
//      out << "\" style=\"fill:none;stroke:red;stroke-width:0.5px\"/>\n";
//    }
  }
}

void writeHeader(std::ofstream &out, o2::mch::contour::BBox<double> &box, int scale)
{
  out << "<svg height=\"" << scale * box.height() + 100 << "\" width=\"" << scale * box.width() + 100 << "\">\n";
}

void writeContour(std::string filename, const o2::mch::contour::Contour<double> &contour, int scale)
{
  std::ofstream out(filename);
  out << "<html><body>\n";
  auto box = getBBox(contour);
  writeHeader(out, box, scale);
  writeContour(out, scale, contour, box);
  out << "</svg>\n";
  out << "</body></html>\n";
}

void writeContours(const std::vector<o2::mch::contour::Contour<double>>& contours, const char *filename, double x,
                   double y)
{
  std::ofstream out(filename);
  int scale = 10;
  out << "<html><body>\n";
  auto env = o2::mch::contour::getEnvelop(contours);
  auto box = getBBox(env);

  x -= box.xmin();
  y -= box.ymin();

  writeHeader(out, box, scale);
  writeContour(out, scale, env, box, "fill:#eeeeee;stroke:black;stroke-width:3px");
  for (auto i = 0; i < contours.size(); ++i) {
    auto &c = contours[i];
//  for (auto &&c: contours) {
    writeContour(out, scale, c, box);
    if (c.contains(x, y)) {
      writeContour(out, scale, c, box, "fill:#aaaaaa;stroke:none");
    }
  }

  out << "<circle cx=\"" << scale * x << "\" cy=\"" << scale * y << "\" r=\"5\"\n"
    "style=\"fill:none;stroke:black;stroke-width:0.5px;\"/>";

  out << "<circle cx=\"" << scale * x << "\" cy=\"" << scale * y << "\" r=\"1\"\n"
    "style=\"fill:none;stroke:red;stroke-width:0.5px;\"/>";

  out << "</svg>\n";
  out << "</body></html>\n";
}


}
}
}
