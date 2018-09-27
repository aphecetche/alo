#include "Digit_generated.h"
#include "MCHMappingInterface/Segmentation.h"
#include <fstream>
#include <iostream>
int readDigits(const char *filename) {

  std::ifstream in(filename);

  int size;

  while (!in.eof()) {
    in.read(reinterpret_cast<char *>(&size), sizeof(int));
    std::cout << "size=" << size << "\n";
    char *buf = new char[size];
    if (in.eof()) {
      continue;
    }

    in.read(buf, size);
    auto digitDE = o2::mch::GetDigitDE(buf);
    int detElemId = digitDE->detElemId();
    std::cout << "DE " << detElemId << "\n";
    auto digitTB = digitDE->digitTimeBlocks()->Get(0);
    auto digitPlane = digitTB->digitPlanes()->Get(0);
    bool isBending = digitPlane->isBending();
    std::cout << ( isBending ? "B":"NB") << "\n";
    auto digits = digitPlane->digits();
    o2::mch::mapping::Segmentation seg{detElemId,isBending};
    for (auto d : *digits) {
        int dsId = seg.padDualSampaId(d->uid());
        int dsCh = seg.padDualSampaChannel(d->uid());
        std::cout << "\t\t" << d->uid() << " [" << dsId << "," << dsCh << "]" << d->adc() << "\n";
    }
    delete[] buf;
  }
}

int main(int argc, char **argv) { return readDigits(argv[1]); }
