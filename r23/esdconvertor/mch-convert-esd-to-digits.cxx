#include "AliESDEvent.h"
#include "AliESDMuonCluster.h"
#include "AliESDMuonPad.h"
#include "AliESDMuonTrack.h"
#include "AliMUONVDigit.h"
#include "Digit_generated.h"
#include "FileCreation.h"
#include "SegmentationPair.h"
#include "TFile.h"
#include "TTree.h"
#include "boost/program_options.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <memory>
#include <vector>

namespace po = boost::program_options;
using namespace o2::mch;

std::vector<AliESDMuonTrack*> getTracks(AliESDEvent& event)
{

  std::vector<AliESDMuonTrack*> tracks;

  auto nTracks = event.GetNumberOfMuonTracks();

  for (auto iTrack = 0; iTrack < nTracks; iTrack++) {
    AliESDMuonTrack* track = event.GetMuonTrack(iTrack);
    if (track->ContainTrackerData()) {
      tracks.push_back(track);
    }
  }
  return tracks;
}

std::vector<AliESDMuonCluster*> getClusters(AliESDEvent& event)
{

  std::vector<AliESDMuonCluster*> clusters;

  auto tracks = getTracks(event);

  for (auto track : tracks) {

    auto nClusters = track->GetNClusters();

    for (auto iCl = 0; iCl < nClusters; iCl++) {
      clusters.push_back(event.FindMuonCluster(track->GetClusterId(iCl)));
    }
  }
  return clusters;
}

std::set<unsigned long> getPadIds(AliESDEvent& event)
{

  std::set<unsigned long> padIds;

  auto clusters = getClusters(event);

  for (auto c : clusters) {
    for (int i = 0; i < c->GetNPads(); ++i) {
      padIds.insert(c->GetPadId(i));
    }
  }

  return padIds;
}

void convertOneDE(const SegmentationPair& seg, int detElemId,
                  AliESDEvent& event, const std::set<unsigned long>& padIds,
                  std::ofstream& out)
{

  std::vector<flatbuffers::Offset<o2::mch::Digit>> bendingDigits;
  std::vector<flatbuffers::Offset<o2::mch::Digit>> nonBendingDigits;
  flatbuffers::FlatBufferBuilder fbb{
    1024
  }; // TODO what default initial size to use ?

  for (auto padid : padIds) {
    int de, manuId, manuChannel, cathode;
    AliMUONVDigit::DecodeUniqueID(padid, de, manuId, manuChannel, cathode);
    if (de != detElemId) {
      continue;
    }
    bool isBending = ((manuId & 1024) == 0);
    auto uid = seg[isBending].findPadByFEE(manuId, manuChannel);
    if (!(seg[isBending].isValid(uid))) {
      std::cout << "got invalid pad !\n";
      std::cout << "[DE" << detElemId << " MANU" << manuId << " CH "
                << manuChannel << "]\n";
    } else {
      auto pad = event.FindMuonPad(padid);
      auto adc = pad->GetADC();
      if (isBending) {
        bendingDigits.push_back(o2::mch::CreateDigit(fbb, uid, adc));
      } else {
        nonBendingDigits.push_back(o2::mch::CreateDigit(fbb, uid, adc));
      }
    }

    // std::cout << "\t\t[DE" << detElemId << " MANU" << manuId << " CH "
    //           << manuChannel << "]";
    // auto pad = event.FindMuonPad(padid);
    // std::cout << pad->GetADC() << "\n";
  }

  if (bendingDigits.empty() && nonBendingDigits.empty()) {
    return;
  }

  flatbuffers::Offset<DigitPlane> b =
    o2::mch::CreateDigitPlaneDirect(fbb, true, &bendingDigits);
  flatbuffers::Offset<DigitPlane> nb =
    o2::mch::CreateDigitPlaneDirect(fbb, false, &nonBendingDigits);

  auto vb = fbb.CreateVector(&b, 1);
  auto vnb = fbb.CreateVector(&nb, 1);

  int timestamp = 0;
  auto tb_b = o2::mch::CreateDigitTimeBlock(fbb, timestamp, vb);
  auto tb_nb = o2::mch::CreateDigitTimeBlock(fbb, timestamp, vnb);

  std::vector<flatbuffers::Offset<DigitTimeBlock>> vtb{ tb_b, tb_nb };

  auto digitDE = o2::mch::CreateDigitDE(fbb, detElemId, fbb.CreateVector(vtb));

  fbb.Finish(digitDE);

  // This must be called after `Finish()`.
  uint8_t* buf = fbb.GetBufferPointer();
  int size = fbb.GetSize(); // Returns the size of the buffer that
                            // `GetBufferPointer()` points to.
  // std::cout << "size=" << size << " b=" << bendingDigits.size()
  //           << " nb=" << nonBendingDigits.size() << "\n";
  out.write(reinterpret_cast<const char*>(&size), sizeof(int));
  out.write(reinterpret_cast<const char*>(buf), sizeof(uint8_t) * size);
}

void convertESD(const std::vector<int>& detElemIds, const char* esdFileName, std::map<int,std::ofstream>& out)
{
  // open the ESD file
  std::unique_ptr<TFile> esdFile{ TFile::Open(esdFileName) };

  if (!esdFile || !esdFile->IsOpen()) {
    Error("convertESD", "opening ESD file %s failed", esdFileName);
    return;
  }

  TTree* tree = static_cast<TTree*>(esdFile->Get("esdTree"));
  if (!tree) {
    Error("convertESD", "no ESD tree found");
    return;
  }

  AliESDEvent event;
  event.ReadFromTree(tree);

  auto nevents = tree->GetEntries();

  std::cout << nevents << " entries to process\n";

  tree->SetBranchStatus("*", false);
  tree->SetBranchStatus("AliESDRun*", true);
  tree->SetBranchStatus("AliESDHeader*", true);
  tree->SetBranchStatus("Muon*", true);

  std::map<int, std::unique_ptr<SegmentationPair>> segmentations;

  for (auto detElemId : detElemIds) {
    segmentations.insert(std::make_pair(
      detElemId, std::make_unique<SegmentationPair>(detElemId)));
  }

  for (auto iEvent = 0; iEvent < nevents; iEvent++) {

    if (tree->GetEvent(iEvent) <= 0) {
      Error("convertESD", "no ESD object found for event %d", iEvent);
      return;
    }

    auto padIds = getPadIds(event);
    if (!padIds.size()) {
      continue;
    }

    for (auto detElemId : detElemIds) {
      auto file = out.find(detElemId);
      convertOneDE(*(segmentations[detElemId]), detElemId, event, padIds,file->second);
    }
  }
}

int main(int argc, char** argv)
{
  po::variables_map vm;
  po::options_description generic("Generic options");
  std::string basename;
  std::vector<int> detElemIds;

  // clang-format off
  generic.add_options()
          ("help", "produce help message")
          ("detelemids", po::value<std::vector<int>>(&detElemIds), "produce digit for those detetcion elements (all if not specified)")
          ("basename", po::value<std::string>(&basename),"basename of output files");

  po::options_description hidden("hidden options");
  hidden.add_options()
          ("input-file", po::value<std::vector<std::string>>(),"input file");
  // clang-format on

  po::options_description cmdline;
  cmdline.add(generic).add(hidden);

  po::positional_options_description p;
  p.add("input-file", -1);

  po::store(
    po::command_line_parser(argc, argv).options(cmdline).positional(p).run(),
    vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << generic << std::endl;
    return 2;
  }
  if (vm.count("input-file") == 0) {
    std::cout << "no input file specified" << std::endl;
    return 1;
  }

  std::vector<std::string> inputfiles{
    vm["input-file"].as<std::vector<std::string>>()
  };

  if (detElemIds.empty()) {
    std::cout << "no detection element specified : using all of them\n";
    o2::mch::mapping::forEachDetectionElement([&detElemIds](int detElemId) { detElemIds.push_back(detElemId); });
  }

  auto out = createDetElemFiles(basename.c_str(), detElemIds);

  for (auto input : inputfiles) {
    convertESD(detElemIds, input.c_str(), out);
  }
  return 0;
}

