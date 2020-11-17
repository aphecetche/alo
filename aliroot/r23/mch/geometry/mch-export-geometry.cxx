#include <boost/program_options.hpp>
#include <iostream>
#include <stdexcept>
#include <TGeoManager.h>
#include <TFile.h>
#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "TClonesArray.h"
#include "AliAlignObjMatrix.h"
#include <sstream>
#include <vector>
#include <string>
#include "TGeoPhysicalNode.h"
#include <rapidjson/document.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

namespace po = boost::program_options;

std::vector<std::string> splitString(const std::string& src, char delim)
{
  std::stringstream ss(src);
  std::string token;
  std::vector<std::string> tokens;

  while (std::getline(ss, token, delim)) {
    if (!token.empty()) {
      tokens.push_back(std::move(token));
    }
  }

  return tokens;
}

TGeoManager* readFromFile(std::string filename)
{
  TFile* f = TFile::Open(filename.c_str());
  if (f->IsZombie()) {
    throw std::runtime_error("can not open " + filename);
  }

  TGeoManager* geo = static_cast<TGeoManager*>(f->Get("ALICE"));
  if (!geo) {
    geo = static_cast<TGeoManager*>(f->Get("FAIRGeom"));
  }
  if (!geo) {
    f->ls();
    throw std::runtime_error("could not find ALICE geometry (using ALICE or FAIRGeom names)");
  }
  return geo;
}

TGeoManager* readFromCCDB(std::string ocdb, int run)
{
  return nullptr;
}

bool align(TGeoManager& geom, std::string ocdb, int run)
{
  auto man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdb.c_str());
  man->SetRun(run);
  auto entry = man->Get("MUON/Align/Data", run);
  if (!entry) {
    throw std::runtime_error("could not get MUON/Align/Data");
  }
  TClonesArray* array = dynamic_cast<TClonesArray*>(entry->GetObject());
  for (int i = 0; i < array->GetEntries(); i++) {
    AliAlignObjMatrix* mat = static_cast<AliAlignObjMatrix*>(array->At(i));

    std::string name = mat->GetSymName();
    auto parts = splitString(name, '/');
    auto id = std::stoi(parts[1].substr(2));
    bool isMCH = id <= 15;
    if (isMCH) {
      std::cout << "Would need to align " << name << " " << id << "\n";
    }
  }
  return false;
}
template <typename WRITER>
void matrix2json(std::string name, const TGeoHMatrix& matrix, WRITER& w)
{

  const Double_t* t = matrix.GetTranslation();
  const Double_t* m = matrix.GetRotationMatrix();
  w.StartObject();
  w.Key("name");
  w.String(name.c_str());
  w.Key("tx");
  w.Double(t[0]);
  w.Key("ty");
  w.Double(t[1]);
  w.Key("tz");
  w.Double(t[2]);
  w.Key("rotation");
  w.StartArray();
  for (auto i = 0; i < 3; i++) {
    w.StartArray();
    w.Double(m[i * 3]);
    w.Double(m[i * 3 + 1]);
    w.Double(m[i * 3 + 2]);
    w.EndArray();
  }
  w.EndArray();
  w.EndObject();
}

bool isMCH(std::string alignableName)
{
  auto parts = splitString(alignableName, '/');
  bool aliroot = parts[0] == "MUON";
  bool o2 = parts[0] == "MCH";
  if (!o2 && !aliroot) {
    return false;
  }
  auto id = std::stoi(parts[1].substr(2));
  return (aliroot && (id <= 15)) || (o2 && (id <= 19));
}

void exportGeom(const TGeoManager& geom)
{
  // export geometry as json document

  rapidjson::OStreamWrapper osw(std::cout);
  rapidjson::Writer<rapidjson::OStreamWrapper> writer(osw);

  writer.StartObject();
  writer.Key("alignables");
  writer.StartArray();

  for (auto i = 0; i < geom.GetNAlignable(); i++) {
    auto ae = geom.GetAlignableEntry(i);
    std::string name = ae->GetName();
    if (!isMCH(ae->GetName())) {
      continue;
    }
    if (ae->GetMatrix()) {
      //ae->GetMatrix()->Print();
      matrix2json(ae->GetName(), *ae->GetMatrix(), writer);
    }
    if (ae->GetMatrixOrig()) {
      std::cout << "WARNING : got also OrigMatrix !!";
      //ae->GetMatrixOrig()->Print();
    }
  }

  writer.EndArray();
  writer.EndObject();
}

int main(int argc, char** argv)
{
  po::variables_map vm;
  po::options_description options;

  // clang-format off
    options.add_options()
     ("help,h","help")
     ("ocdb",po::value<std::string>(),"OCDB path for Run2 geometry")
     ("geom",po::value<std::string>(),"geometry.root file")
     ("run",po::value<int>(),"run number");
  // clang-format on

  po::options_description cmdline;
  cmdline.add(options);

  po::store(po::command_line_parser(argc, argv).options(cmdline).run(), vm);

  if (vm.count("help")) {
    std::cout << "This program extracts MCH geometry\n";
  }

  TGeoManager* geom{ nullptr };

  if (vm.count("geom")) {
    geom = readFromFile(vm["geom"].as<std::string>());
  }

  if (vm.count("run") && vm.count("ocdb")) {
    auto run = vm["run"].as<int>();
    auto ocdb = vm["ocdb"].as<std::string>();
    if (!geom) {
      geom = readFromCCDB(ocdb, run);
    }
    align(*geom, ocdb, run);
  }

  exportGeom(*geom);

  return 0;
}
