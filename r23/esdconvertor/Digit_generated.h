// automatically generated by the FlatBuffers compiler, do not modify


#ifndef FLATBUFFERS_GENERATED_DIGIT_O2_MCH_H_
#define FLATBUFFERS_GENERATED_DIGIT_O2_MCH_H_

#include "flatbuffers/flatbuffers.h"

namespace o2 {
namespace mch {

struct Digit;

struct DigitEvent;

struct Digit FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
  enum {
    VT_UID = 4,
    VT_ADC = 6
  };
  int32_t uid() const {
    return GetField<int32_t>(VT_UID, 0);
  }
  int32_t adc() const {
    return GetField<int32_t>(VT_ADC, 0);
  }
  bool Verify(flatbuffers::Verifier &verifier) const {
    return VerifyTableStart(verifier) &&
           VerifyField<int32_t>(verifier, VT_UID) &&
           VerifyField<int32_t>(verifier, VT_ADC) &&
           verifier.EndTable();
  }
};

struct DigitBuilder {
  flatbuffers::FlatBufferBuilder &fbb_;
  flatbuffers::uoffset_t start_;
  void add_uid(int32_t uid) {
    fbb_.AddElement<int32_t>(Digit::VT_UID, uid, 0);
  }
  void add_adc(int32_t adc) {
    fbb_.AddElement<int32_t>(Digit::VT_ADC, adc, 0);
  }
  explicit DigitBuilder(flatbuffers::FlatBufferBuilder &_fbb)
        : fbb_(_fbb) {
    start_ = fbb_.StartTable();
  }
  DigitBuilder &operator=(const DigitBuilder &);
  flatbuffers::Offset<Digit> Finish() {
    const auto end = fbb_.EndTable(start_);
    auto o = flatbuffers::Offset<Digit>(end);
    return o;
  }
};

inline flatbuffers::Offset<Digit> CreateDigit(
    flatbuffers::FlatBufferBuilder &_fbb,
    int32_t uid = 0,
    int32_t adc = 0) {
  DigitBuilder builder_(_fbb);
  builder_.add_adc(adc);
  builder_.add_uid(uid);
  return builder_.Finish();
}

struct DigitEvent FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
  enum {
    VT_BENDINGDIGITS = 4,
    VT_NONBENDINGDIGITS = 6
  };
  const flatbuffers::Vector<flatbuffers::Offset<Digit>> *bendingDigits() const {
    return GetPointer<const flatbuffers::Vector<flatbuffers::Offset<Digit>> *>(VT_BENDINGDIGITS);
  }
  const flatbuffers::Vector<flatbuffers::Offset<Digit>> *nonBendingDigits() const {
    return GetPointer<const flatbuffers::Vector<flatbuffers::Offset<Digit>> *>(VT_NONBENDINGDIGITS);
  }
  bool Verify(flatbuffers::Verifier &verifier) const {
    return VerifyTableStart(verifier) &&
           VerifyOffset(verifier, VT_BENDINGDIGITS) &&
           verifier.Verify(bendingDigits()) &&
           verifier.VerifyVectorOfTables(bendingDigits()) &&
           VerifyOffset(verifier, VT_NONBENDINGDIGITS) &&
           verifier.Verify(nonBendingDigits()) &&
           verifier.VerifyVectorOfTables(nonBendingDigits()) &&
           verifier.EndTable();
  }
};

struct DigitEventBuilder {
  flatbuffers::FlatBufferBuilder &fbb_;
  flatbuffers::uoffset_t start_;
  void add_bendingDigits(flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<Digit>>> bendingDigits) {
    fbb_.AddOffset(DigitEvent::VT_BENDINGDIGITS, bendingDigits);
  }
  void add_nonBendingDigits(flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<Digit>>> nonBendingDigits) {
    fbb_.AddOffset(DigitEvent::VT_NONBENDINGDIGITS, nonBendingDigits);
  }
  explicit DigitEventBuilder(flatbuffers::FlatBufferBuilder &_fbb)
        : fbb_(_fbb) {
    start_ = fbb_.StartTable();
  }
  DigitEventBuilder &operator=(const DigitEventBuilder &);
  flatbuffers::Offset<DigitEvent> Finish() {
    const auto end = fbb_.EndTable(start_);
    auto o = flatbuffers::Offset<DigitEvent>(end);
    return o;
  }
};

inline flatbuffers::Offset<DigitEvent> CreateDigitEvent(
    flatbuffers::FlatBufferBuilder &_fbb,
    flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<Digit>>> bendingDigits = 0,
    flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<Digit>>> nonBendingDigits = 0) {
  DigitEventBuilder builder_(_fbb);
  builder_.add_nonBendingDigits(nonBendingDigits);
  builder_.add_bendingDigits(bendingDigits);
  return builder_.Finish();
}

inline flatbuffers::Offset<DigitEvent> CreateDigitEventDirect(
    flatbuffers::FlatBufferBuilder &_fbb,
    const std::vector<flatbuffers::Offset<Digit>> *bendingDigits = nullptr,
    const std::vector<flatbuffers::Offset<Digit>> *nonBendingDigits = nullptr) {
  return o2::mch::CreateDigitEvent(
      _fbb,
      bendingDigits ? _fbb.CreateVector<flatbuffers::Offset<Digit>>(*bendingDigits) : 0,
      nonBendingDigits ? _fbb.CreateVector<flatbuffers::Offset<Digit>>(*nonBendingDigits) : 0);
}

inline const o2::mch::DigitEvent *GetDigitEvent(const void *buf) {
  return flatbuffers::GetRoot<o2::mch::DigitEvent>(buf);
}

inline bool VerifyDigitEventBuffer(
    flatbuffers::Verifier &verifier) {
  return verifier.VerifyBuffer<o2::mch::DigitEvent>(nullptr);
}

inline void FinishDigitEventBuffer(
    flatbuffers::FlatBufferBuilder &fbb,
    flatbuffers::Offset<o2::mch::DigitEvent> root) {
  fbb.Finish(root);
}

}  // namespace mch
}  // namespace o2

#endif  // FLATBUFFERS_GENERATED_DIGIT_O2_MCH_H_
