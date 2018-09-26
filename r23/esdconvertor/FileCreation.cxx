#include <fstream>
#include <map>
#include <vector>
#include <string>

std::map<int, std::ofstream> createDetElemFiles(const char *baseName, const std::vector<int>& detElemIds) {

  std::map<int, std::ofstream> outputFiles;

  for (auto detElemId: detElemIds) {
    outputFiles[detElemId].open(baseName + std::to_string(detElemId) + ".dat");
  }

  return outputFiles;
}
