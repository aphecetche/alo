#include "AliHLTMUONMapping.h"

#include "AliCDBManager.h"
#include "AliMpCDB.h"
#include "AliMpDEIterator.h"
#include "AliMpSegmentation.h"
#include "AliMpVPadIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONVDigit.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

UShort_t GetColumn ( Int_t localBoardId );

std::vector<AliHLTMUONMapping::mpDE> CreateMapping( Int_t firstChamber, Int_t lastChamber );

void FindNeighbours(AliHLTMUONMapping::mpDE &de, UChar_t iPlane);

int WriteMapping(std::vector<AliHLTMUONMapping::mpDE> detectionElements, const char* mapfile);

//______________________________________________________________________________
UShort_t GetColumn ( Int_t localBoardId )
{
  /// Get column from local board id
  Int_t iboard = (localBoardId-1)%117;
  Int_t colChange[6] = {16,38,60,76,92,108};
  UShort_t icol = 6;
  for ( UShort_t ichange=0; ichange<6; ichange++ ) {
    if ( iboard < colChange[ichange] ) {
      icol = ichange;
      break;
    }
  }
  return icol;
}

//______________________________________________________________________________
std::vector<AliHLTMUONMapping::mpDE> CreateMapping( Int_t firstChamber, Int_t lastChamber )
{
  /// fill the mpDE and mpPad structures once for all

  // loop over DEs

  cout << "Creating mapping..." << endl;

  std::vector<AliHLTMUONMapping::mpDE> detectionElements;

  UInt_t padId = 0;
  for (Int_t iCh = firstChamber; iCh <= lastChamber; ++iCh) {
    AliMpDEIterator deIt;
    deIt.First(iCh);

    while (!deIt.IsDone()) {

      Int_t deId = deIt.CurrentDEId();

      cout << "DE " << deId << "..." << flush;

      AliHLTMUONMapping::mpDE de;
      de.id = deId;

      const AliMpVSegmentation* seg[2] =
      { AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath0),
        AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath1)
      };

      UChar_t iPlane[2];
      if (seg[0]->PlaneType() == AliMp::kNonBendingPlane) {
        iPlane[0] = 1;
        iPlane[1] = 0;
      } else {
        iPlane[0] = 0;
        iPlane[1] = 1;
      }

      de.nPads[iPlane[0]] = seg[0]->NofPads();
      de.nPads[iPlane[1]] = seg[1]->NofPads();
      de.pads = new AliHLTMUONMapping::mpPad[de.nPads[0]+de.nPads[1]];

      // de.ownDigits.reserve((de.nPads[0]/10+de.nPads[1]/10)); // 10% occupancy
      // de.nOwnDigits = 0;
      // de.digits.reserve((de.nPads[0]/10+de.nPads[1]/10)); // 10% occupancy
      // de.nOrderedPads[0] = 0;
      // de.orderedPads[0].reserve((de.nPads[0]/10+de.nPads[1]/10)); // 10% occupancy
      // de.nOrderedPads[1] = 0;
      // de.orderedPads[1].reserve((de.nPads[0]/10+de.nPads[1]/10)); // 10% occupancy

      // loop over cathods
      for (UChar_t iCath = 0; iCath < 2; ++iCath) {

        de.iCath[iPlane[iCath]] = iCath;

        // de.nFiredPads[iPlane[iCath]] = 0;
        // de.firedPads[iPlane[iCath]].reserve(de.nPads[iPlane[iCath]]/10); // 10% occupancy

        // loop over pads to associate an index to their Id
        UShort_t iPad = (iPlane[iCath] == 0) ? 0 : de.nPads[0];
        AliMpVPadIterator* padIt = seg[iCath]->CreateIterator();
        padIt->First();
        while (!padIt->IsDone()) {

          AliMpPad pad = padIt->CurrentItem();

          padId = AliMUONVDigit::BuildUniqueID(deId, pad.GetManuId(), pad.GetManuChannel(), iCath);
          de.padIndices[iPlane[iCath]].Add(padId, iPad+1);

          de.pads[iPad].iDigit = 0;
          if ( deId >= 1100 ) de.pads[iPad].iDigit = GetColumn(pad.GetManuId());
          de.pads[iPad].nNeighbours = 0;
          de.pads[iPad].area[0][0] = pad.GetPositionX() - pad.GetDimensionX();
          de.pads[iPad].area[0][1] = pad.GetPositionX() + pad.GetDimensionX();
          de.pads[iPad].area[1][0] = pad.GetPositionY() - pad.GetDimensionY();
          de.pads[iPad].area[1][1] = pad.GetPositionY() + pad.GetDimensionY();
          de.pads[iPad].useMe = kFALSE;

          ++iPad;
          padIt->Next();

        }

        FindNeighbours(de, iPlane[iCath]);

        delete padIt;

      }

      deIt.Next();

      detectionElements.push_back(de);

    }

  }

  return detectionElements;
}

//______________________________________________________________________________
void FindNeighbours(AliHLTMUONMapping::mpDE &de, UChar_t iPlane)
{
  /// fill the mpPad neighbours structures of given DE/plane

  UShort_t firstPad = (iPlane == 0) ? 0 : de.nPads[0];
  UShort_t lastPad = firstPad + de.nPads[iPlane];

  // loop over pads
  for (UShort_t iPad1 = firstPad; iPad1 < lastPad; ++iPad1) {

    // loop over next pads to find the neighbours
    for (UShort_t iPad2 = iPad1+1; iPad2 < lastPad; ++iPad2) {

      // For the tracking, the pad neighbours of x are as following:
      //  --- --- ---
      // |   |   |   |
      //  --- --- ---
      // |   | x |   |
      //  --- --- ---
      // |   |   |   |
      //  --- --- ---
      // For the trigger, which has much larger strip, we define as neighbours:
      //      ---
      //     |   |
      //  --- --- ---
      // |   | x |   |
      //  --- --- ---
      //     |   |
      //      ---

      Bool_t hasOverlap = ( de.id >= 1100 ) ? AliHLTMUONMapping::AreOverlappingExcludeCorners(de.pads[iPad1].area, de.pads[iPad2].area) : AliHLTMUONMapping::AreOverlapping(de.pads[iPad1].area, de.pads[iPad2].area, 1.e-4);
      if ( hasOverlap ) {

        de.pads[iPad1].neighbours[de.pads[iPad1].nNeighbours] = iPad2;
        ++de.pads[iPad1].nNeighbours;
        de.pads[iPad2].neighbours[de.pads[iPad2].nNeighbours] = iPad1;
        ++de.pads[iPad2].nNeighbours;

      }

    }

  }

}

//______________________________________________________________________________
int WriteMapping(std::vector<AliHLTMUONMapping::mpDE> detectionElements, const char* mapfile)
{
  cout << "Writing mapping to " << mapfile << endl;

  ofstream out(mapfile,ios::out | ios::binary);
  if (!out.is_open()) {
    return 1;
  }

  Int_t numberOfDetectionElements = detectionElements.size();

  out.write((char*)&numberOfDetectionElements,sizeof(Int_t));

  for ( AliHLTMUONMapping::mpDE &de : detectionElements ) {

    out.write((char*)&de.id,sizeof(Int_t));
    out.write((char*)&de.iCath[0],sizeof(UChar_t)*2);
    out.write((char*)&de.nPads[0],sizeof(UShort_t)*2);

    Int_t n = de.nPads[0] + de.nPads[1];
    for ( int ip = 0; ip < n; ++ip ) {
      out.write((char*)&(de.pads[ip]),sizeof(struct AliHLTMUONMapping::mpPad));
    }


    Long64_t key, value;
    Int_t mapsize = de.padIndices[0].GetSize() + de.padIndices[1].GetSize();
    assert( n == mapsize);
    mapsize *= 2;
    Long64_t* themap = new Long64_t[mapsize];
    int ix(0);

    for ( int cathode = 0; cathode < 2; ++cathode) {
      TExMapIter mit(&de.padIndices[cathode]);

      while (mit.Next(key,value)) {
          themap[ix++] = key;
          themap[ix++] = value;
      }
    }

    assert(ix==mapsize);

    out.write((char*)&themap[0],sizeof(Long64_t)*mapsize);

    delete[] themap;
  }

  return 0;
}

//______________________________________________________________________________
int main(int argc, char** argv) {

  const char* cdbPath = NULL;
  Int_t run = -1;
  string binmapfile="binmapfile.dat";
  Int_t chamberRange[2] = {0,AliMUONConstants::NCh()-1};

  for (int i = 0; i < argc; i++) {

    // get the binary mapping file name
    if  (strcmp(argv[i],"-binmapfile")==0) {

      if (!binmapfile.empty()) {
        cout << "binmapfile already specified. Will replace previous value" << endl;
      }

      if  (argc <= i+1 ) {
        cerr << "binmapfile was not specified" << endl;
        return 1;
      }
      binmapfile = argv[i+1];
      ++i;
      continue;
    }

    // get the CDB path
    if (strcmp(argv[i], "-cdbpath") == 0) {

      if (cdbPath != NULL) {
        cout << "CDB path was already specified. Will replace previous value given by -cdbpath." << endl;
      }

      if (argc <= i+1) {
        cerr << "The CDB path was not specified." << endl;
        return 2;
      }

      cdbPath = argv[i+1];
      ++i;
      continue;
    }

    // get the run number
    if (strcmp(argv[i], "-run") == 0) {

      if (run != -1) {
        cout << "Run number was already specified. Will replace previous value given by -run." << endl;
      }

      if (argc <= i+1) {
        cerr << "The run number was not specified." << endl;
        return 3;
      }

      char* cpErr = NULL;
      run = Int_t( strtol(argv[i+1], &cpErr, 0) );
      if (cpErr == NULL or *cpErr != '\0' or run < 0) {
        cerr << "Cannot convert "<< argv[i+1] << " to a valid run number. Expected a positive integer value." << endl;
        return 4;
      }

      cout << "run=" << run << endl;

      ++i;
      continue;
    }

    // get the chamber range
    for ( Int_t ival=0; ival<2; ival++ ) {
      const char* argName = ( ival==0 ) ? "-first" : "-last";
      if (strcmp(argv[i], argName) == 0) {

        if (argc <= i+1) {
          cerr << "The run number was not specified." << endl;
          return 3;
        }

        char* cpErr = NULL;
        chamberRange[ival] = Int_t( strtol(argv[i+1], &cpErr, 0) );
        if (cpErr == NULL or *cpErr != '\0' or run < 0) {
          cerr << "Cannot convert "<< argv[i+1] << " to a valid run number. Expected a positive integer value." << endl;
          return 4;
        }

        cout << argName << chamberRange[ival] << endl;

        ++i;
        continue;
      }
    }
  }


  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if (cdbManager == NULL) {
    cerr << "CDB manager instance does not exist." << endl;
    return 4;
  }

  // Set the CDB path
  if (cdbPath == NULL) {
    cerr << "The CDB path was not set." << endl;
    return 5;
  } else {
    cdbManager->SetDefaultStorage(cdbPath);
    if (!cdbManager->IsDefaultStorageSet()) {
      cerr << "Invalid CDB path: " << cdbPath << endl;
      return 6;
    }
  }

  // Set the run number
  if (run == -1) {
    cerr << "The run number was not set." << endl;
    return 7;
  } else {
    cdbManager->SetRun(run);
  }

  // Check the chamber range
  if ( chamberRange[0] > chamberRange[1] ) {
    cerr << "First chamber (" << chamberRange[0] << ") must be smaller than last chamber (" << chamberRange[1] << ")" << endl;
  }
  else {
    for ( Int_t ival=0; ival>2; ival++ ) {
      if ( chamberRange[ival] >= 0 && chamberRange[ival] < AliMUONConstants::NCh() ) continue;
      cerr << "Chamber number must be in the range 0 - 13. MCH: 0-9; MTR: 10-13" << endl;
      return 8;
    }
  }

  // Load the MUON mapping from CDB and create the internal mapping
  // MUON mapping also needed for raw data decoding
  if (!AliMUONCDB::LoadMapping()) return 8;

  std::vector<AliHLTMUONMapping::mpDE> detectionElements = CreateMapping(chamberRange[0],chamberRange[1]);

  int rv = WriteMapping(detectionElements, binmapfile.c_str());

  return rv;
}
