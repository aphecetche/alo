#include "DumpDigits.h"
#include "Run3_generated.h"
#include "SegmentationPair.h"
#include <fstream>
#include <iostream>
#include <map>

void dumpDigitPlane(const o2::mch::DigitPlane* digitPlane, const SegmentationPair& segpair, unsigned long& npadids, unsigned long& nduppadids)
{
  bool isBending = digitPlane->isBending();
  //   std::cout << (isBending ? "B" : "NB") << "\n";
  auto digits = digitPlane->digits();
  auto& seg = segpair[isBending];
  std::set<int> uniq;
  for (auto d : *digits) {
    ++npadids;
    int dsId = seg.padDualSampaId(d->uid());
    int dsCh = seg.padDualSampaChannel(d->uid());
    //     std::cout << "\t\t" << d->uid() << " [" << dsId << "," << dsCh << "]" << d->adc() << "\n";
    uniq.insert(d->uid());
  }
  nduppadids += digits->size() - uniq.size();

  if (uniq.size() != digits->size()) {
    //std::cout << "Got duplicates " << uniq.size() << " " << digits->Length() << "\n";
  }
}

void dumpDigits(const std::string& filename)
{
  std::ifstream in(filename);

  int size;
  int nevents{ 0 };
  unsigned long npadids{ 0 }, nduppadids{ 0 };

  while (!in.eof()) {
    in.read(reinterpret_cast<char*>(&size), sizeof(int));
    //     std::cout << "size=" << size << "\n";
    char* buf = new char[size];
    if (in.eof()) {
      continue;
    }

    nevents++;

    in.read(buf, size);

    auto digitDE = flatbuffers::GetRoot<o2::mch::DigitDE>(buf);

    int detElemId = digitDE->detElemId();

    int nblocks = digitDE->digitTimeBlocks()->size();

    //   std::cout << "DE " << detElemId << " " << nblocks << " blocks\n";

    for (auto i = 0; i < nblocks; i++) {
      auto digitTB = digitDE->digitTimeBlocks()->Get(i);
      for (auto j = 0; j < digitTB->digitPlanes()->size(); ++j) {
        auto digitPlane = digitTB->digitPlanes()->Get(j);
        auto& segpair = getSegmentationPair(detElemId);
        dumpDigitPlane(digitPlane, segpair, npadids, nduppadids);
      }
    }
    delete[] buf;
  }

  std::cout << "nevents=" << nevents << std::endl;
  std::cout << "# of dup padids=" << nduppadids << std::endl;
  std::cout << "# of padids=" << npadids << std::endl;
}
