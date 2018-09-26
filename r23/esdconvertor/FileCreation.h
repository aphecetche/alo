#ifndef FILECREATION_H
#define FILECREATION_H

#include <map>
#include <fstream>
#include <vector>

std::map<int, std::ofstream> createDetElemFiles(const char *baseName, const std::vector<int>& detElemIds);

#endif
