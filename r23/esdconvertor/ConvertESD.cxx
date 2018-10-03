#include "ConvertESD.h"

#include "AliESDEvent.h"
#include "AliESDMuonCluster.h"
#include "AliESDMuonPad.h"
#include "AliESDMuonTrack.h"
#include "Digit_generated.h"
#include "GetESDClusters.h"
#include "SegmentationPair.h"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <memory>
#include <vector>

using namespace o2::mch;

SegmentationMap getSegmentations(const std::vector<int>& detElemIds)
{
  SegmentationMap segmentations;
  for (auto detElemId : detElemIds) {
    segmentations.insert(std::make_pair(
      detElemId, std::make_unique<SegmentationPair>(detElemId)));
  }
  return segmentations;
}

void writeBuffer(flatbuffers::FlatBufferBuilder& fbb, std::ofstream& out)
{
  uint8_t* buf = fbb.GetBufferPointer();
  int size = fbb.GetSize(); // Returns the size of the buffer that
                            // `GetBufferPointer()` points to.
  if (size) {
    out.write(reinterpret_cast<const char*>(&size), sizeof(int));
    out.write(reinterpret_cast<const char*>(buf), sizeof(uint8_t) * size);
  }
}

void convertESD(const SegmentationMap& segmentations,
                const char* esdFileName,
                OutputFileMap& clusterFiles,
                OutputFileMap& digitFiles)
{
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

  flatbuffers::FlatBufferBuilder fbb{ 1024 };

  for (auto iEvent = 0; iEvent < nevents; iEvent++) {

    if (tree->GetEvent(iEvent) <= 0) {
      Error("convertESD", "no ESD object found for event %d", iEvent);
      return;
    }

    for (auto& seg : segmentations) {
      auto detElemId = seg.first;
      auto clusterFile = clusterFiles.find(detElemId);
      auto clusters = getClusters(event, detElemId);
      fbb.Clear();
      convertClusters(*(seg.second), clusters, fbb);
      writeBuffer(fbb, clusterFile->second);
      auto digitFile = digitFiles.find(detElemId);
      if (digitFile->second.is_open() > 0) {
        fbb.Clear();
        convertDigits(*(seg.second), event, clusters, fbb);
        writeBuffer(fbb, digitFile->second);
      }
    }
  }
}

std::vector<unsigned long> getPadIds(std::vector<AliESDMuonCluster*>& clusters)
{
  std::vector<unsigned long> padIds;

  for (auto c : clusters) {
    for (int i = 0; i < c->GetNPads(); ++i) {
      padIds.emplace_back(c->GetPadId(i));
    }
  }

  return padIds;
}

void convertClusters(const SegmentationPair& seg,
                     std::vector<AliESDMuonCluster*>& clusters, flatbuffers::FlatBufferBuilder& fbb)
{
}

void convertDigits(const SegmentationPair& seg,
                   AliESDEvent& event,
                   std::vector<AliESDMuonCluster*>& clusters, flatbuffers::FlatBufferBuilder& fbb)
{

  std::vector<flatbuffers::Offset<o2::mch::Digit>> bendingDigits;
  std::vector<flatbuffers::Offset<o2::mch::Digit>> nonBendingDigits;

  auto padIds = getPadIds(clusters);

  int detElemId{ 0 };

  for (auto padid : padIds) {
    detElemId = padid & 0xFFF;
    int manuId = (padid & 0xFFF000) >> 12;
    int manuChannel = (padid & 0x3F000000) >> 24;
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

    std::cout << "\t\t[DE" << detElemId << " MANU" << manuId << " CH "
              << manuChannel << "]";
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
}
