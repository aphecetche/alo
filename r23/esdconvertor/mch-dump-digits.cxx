#include "Digit_generated.h"
#include "SegmentationPair.h"
#include <fstream>
#include <iostream>
#include <map>

const SegmentationPair& getSegmentationPair(int detElemId)
{
  static std::map<int, std::unique_ptr<SegmentationPair>> segpairs;

  auto it = segpairs.find(detElemId);
  if (it != segpairs.end()) {
    return *(it->second);
  } else {
    segpairs.insert(std::make_pair(detElemId, std::make_unique<SegmentationPair>(detElemId)));
    return getSegmentationPair(detElemId);
  }
}

void dumpDigitPlane(const o2::mch::DigitPlane* digitPlane, const SegmentationPair& segpair)
{
  bool isBending = digitPlane->isBending();
//   std::cout << (isBending ? "B" : "NB") << "\n";
  auto digits = digitPlane->digits();
  auto& seg = segpair[isBending];
  std::set<int> uniq;
  for (auto d : *digits) {
    int dsId = seg.padDualSampaId(d->uid());
    int dsCh = seg.padDualSampaChannel(d->uid());
//     std::cout << "\t\t" << d->uid() << " [" << dsId << "," << dsCh << "]" << d->adc() << "\n";
    uniq.insert(d->uid());
  }
  if (uniq.size() != digits->Length()) {
    std::cout << "Got duplicates " << uniq.size() << " " << digits->Length() << "\n";
  }
}

int readDigits(const char* filename)
{
  std::ifstream in(filename);

  int size;

  while (!in.eof()) {
    in.read(reinterpret_cast<char*>(&size), sizeof(int));
     std::cout << "size=" << size << "\n";
    char* buf = new char[size];
    if (in.eof()) {
      continue;
    }

    in.read(buf, size);

    auto digitDE = o2::mch::GetDigitDE(buf);

    int detElemId = digitDE->detElemId();

    int nblocks = digitDE->digitTimeBlocks()->Length();

     std::cout << "DE " << detElemId << " " << nblocks << " blocks\n";

    for (auto i = 0; i < nblocks; i++) {
      auto digitTB = digitDE->digitTimeBlocks()->Get(i);
      for (auto j = 0; j < digitTB->digitPlanes()->Length(); ++j) {
        auto digitPlane = digitTB->digitPlanes()->Get(j);
        auto& segpair = getSegmentationPair(detElemId);
        dumpDigitPlane(digitPlane, segpair);
      }
    }
    delete[] buf;
  }
}

int main(int argc, char** argv) { return readDigits(argv[1]); }
