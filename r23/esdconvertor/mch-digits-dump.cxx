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
    in.read(buf, size);
    auto digitEvent = o2::mch::GetDigitEvent(buf);
    auto bendingDigits = digitEvent->bendingDigits();
    for (auto d : *bendingDigits) {
        std::cout << "\t\t" << d->uid() << " " << d->adc() << "\n";
    }
    delete[] buf;
  }
}

int main(int argc, char **argv) { return readDigits(argv[1]); }
