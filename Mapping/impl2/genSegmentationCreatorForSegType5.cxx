//
// This file has been generated. Do not modify it by hand or your changes might be lost.
//
#include "segmentationCreator.h"

namespace o2
{
namespace mch
{
namespace mapping
{
namespace impl2
{
Segmentation* createSegType5(bool isBendingPlane)
{
  if (isBendingPlane) {
    return new Segmentation{
      5,
      true,
      /* PG */
      { { 4, 2, 0, 20, -20 },   { 5, 8, 0, 25, -20 },   { 6, 4, 0, 35, -20 },   { 7, 9, 0, 40, -20 },
        { 8, 3, 0, 50, -20 },   { 12, 2, 0, -20, -20 }, { 13, 8, 0, -15, -20 }, { 14, 4, 0, -5, -20 },
        { 15, 9, 0, 0, -20 },   { 16, 3, 0, 10, -20 },  { 20, 2, 0, -60, -20 }, { 21, 8, 0, -55, -20 },
        { 22, 4, 0, -45, -20 }, { 23, 9, 0, -40, -20 }, { 24, 3, 0, -30, -20 }, { 101, 1, 0, 50, 0 },
        { 102, 7, 0, 40, 0 },   { 103, 5, 0, 35, 4 },   { 104, 6, 0, 25, 0 },   { 105, 0, 0, 20, 0 },
        { 110, 1, 0, 10, 0 },   { 111, 7, 0, 0, 0 },    { 112, 5, 0, -5, 4 },   { 113, 6, 0, -15, 0 },
        { 114, 0, 0, -20, 0 },  { 119, 1, 0, -30, 0 },  { 120, 7, 0, -40, 0 },  { 121, 5, 0, -45, 4 },
        { 122, 6, 0, -55, 0 },  { 123, 0, 0, -60, 0 } },
      /* PGT */
      { /* L5 */ { 2, 40, { 55, 56, 54, 57, 53, 58, 52, 59, 51, 60, 50, 61, 49, 62, 48, 63, 31, 32, 30, 33,
                            29, 34, 28, 35, 27, 36, 26, 37, 25, 38, 24, 39, 23, 40, 22, 41, 21, 42, 20, 43,
                            19, 44, 18, 45, 17, 46, 16, 47, 15, -1, 14, -1, 13, -1, 12, -1, 11, -1, 10, -1,
                            9,  -1, 8,  -1, 7,  -1, 6,  -1, 5,  -1, 4,  -1, 3,  -1, 2,  -1, 1,  -1, 0,  -1 } },
        /* L6 */ { 2, 40, { 23, 24, 22, 25, 21, 26, 20, 27, 19, 28, 18, 29, 17, 30, 16, 31, 15, 48, 14, 49,
                            13, 50, 12, 51, 11, 52, 10, 53, 9,  54, 8,  55, 7,  56, 6,  57, 5,  58, 4,  59,
                            3,  60, 2,  61, 1,  62, 0,  63, -1, 32, -1, 33, -1, 34, -1, 35, -1, 36, -1, 37,
                            -1, 38, -1, 39, -1, 40, -1, 41, -1, 42, -1, 43, -1, 44, -1, 45, -1, 46, -1, 47 } },
        /* L7 */ { 2, 40, { 47, -1, 46, -1, 45, -1, 44, -1, 43, -1, 42, -1, 41, -1, 40, -1, 39, -1, 38, -1,
                            37, -1, 36, -1, 35, -1, 34, -1, 33, -1, 32, -1, 63, 0,  62, 1,  61, 2,  60, 3,
                            59, 4,  58, 5,  57, 6,  56, 7,  55, 8,  54, 9,  53, 10, 52, 11, 51, 12, 50, 13,
                            49, 14, 48, 15, 31, 16, 30, 17, 29, 18, 28, 19, 27, 20, 26, 21, 25, 22, 24, 23 } },
        /* L8 */ { 2, 40, { -1, 0,  -1, 1,  -1, 2,  -1, 3,  -1, 4,  -1, 5,  -1, 6,  -1, 7,  -1, 8,  -1, 9,
                            -1, 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15, 47, 16, 46, 17, 45, 18, 44, 19,
                            43, 20, 42, 21, 41, 22, 40, 23, 39, 24, 38, 25, 37, 26, 36, 27, 35, 28, 34, 29,
                            33, 30, 32, 31, 63, 48, 62, 49, 61, 50, 60, 51, 59, 52, 58, 53, 57, 54, 56, 55 } },
        /* O10 */ { 2, 32, { 48, 31, 49, 30, 50, 29, 51, 28, 52, 27, 53, 26, 54, 25, 55, 24, 56, 23, 57, 22, 58, 21,
                             59, 20, 60, 19, 61, 18, 62, 17, 63, 16, 32, 15, 33, 14, 34, 13, 35, 12, 36, 11, 37, 10,
                             38, 9,  39, 8,  40, 7,  41, 6,  42, 5,  43, 4,  44, 3,  45, 2,  46, 1,  47, 0 } },
        /* O9 */ { 2, 32, { 0,  47, 1,  46, 2,  45, 3,  44, 4,  43, 5,  42, 6,  41, 7,  40, 8,  39, 9,  38, 10, 37,
                            11, 36, 12, 35, 13, 34, 14, 33, 15, 32, 16, 63, 17, 62, 18, 61, 19, 60, 20, 59, 21, 58,
                            22, 57, 23, 56, 24, 55, 25, 54, 26, 53, 27, 52, 28, 51, 29, 50, 30, 49, 31, 48 } },
        /* Z1 */ { 3, 40, { -1, 39, 40, -1, 38, 41, -1, 37, 42, -1, 36, 43, -1, 35, 44, -1, 34, 45, -1, 33,
                            46, -1, 32, 47, -1, 63, -1, -1, 62, -1, -1, 61, -1, -1, 60, -1, -1, 59, -1, -1,
                            58, -1, -1, 57, -1, -1, 56, -1, -1, 55, -1, -1, 54, -1, -1, 53, -1, -1, 52, -1,
                            -1, 51, -1, -1, 50, -1, -1, 49, -1, -1, 48, -1, 0,  31, -1, 1,  30, -1, 2,  29,
                            -1, 3,  28, -1, 4,  27, -1, 5,  26, -1, 6,  25, -1, 7,  24, -1, 8,  23, -1, 9,
                            22, -1, 10, 21, -1, 11, 20, -1, 12, 19, -1, 13, 18, -1, 14, 17, -1, 15, 16, -1 } },
        /* Z2 */ { 3, 40, { 7,  8,  -1, 6,  9,  -1, 5,  10, -1, 4,  11, -1, 3,  12, -1, 2,  13, -1, 1,  14,
                            -1, 0,  15, -1, -1, 16, -1, -1, 17, -1, -1, 18, -1, -1, 19, -1, -1, 20, -1, -1,
                            21, -1, -1, 22, -1, -1, 23, -1, -1, 24, -1, -1, 25, -1, -1, 26, -1, -1, 27, -1,
                            -1, 28, -1, -1, 29, -1, -1, 30, -1, -1, 31, -1, -1, 48, 47, -1, 49, 46, -1, 50,
                            45, -1, 51, 44, -1, 52, 43, -1, 53, 42, -1, 54, 41, -1, 55, 40, -1, 56, 39, -1,
                            57, 38, -1, 58, 37, -1, 59, 36, -1, 60, 35, -1, 61, 34, -1, 62, 33, -1, 63, 32 } },
        /* Z3 */ { 3, 40, { 32, 63, -1, 33, 62, -1, 34, 61, -1, 35, 60, -1, 36, 59, -1, 37, 58, -1, 38, 57,
                            -1, 39, 56, -1, 40, 55, -1, 41, 54, -1, 42, 53, -1, 43, 52, -1, 44, 51, -1, 45,
                            50, -1, 46, 49, -1, 47, 48, -1, -1, 31, -1, -1, 30, -1, -1, 29, -1, -1, 28, -1,
                            -1, 27, -1, -1, 26, -1, -1, 25, -1, -1, 24, -1, -1, 23, -1, -1, 22, -1, -1, 21,
                            -1, -1, 20, -1, -1, 19, -1, -1, 18, -1, -1, 17, -1, -1, 16, -1, -1, 15, 0,  -1,
                            14, 1,  -1, 13, 2,  -1, 12, 3,  -1, 11, 4,  -1, 10, 5,  -1, 9,  6,  -1, 8,  7 } },
        /* Z4 */ { 3, 40, { -1, 16, 15, -1, 17, 14, -1, 18, 13, -1, 19, 12, -1, 20, 11, -1, 21, 10, -1, 22,
                            9,  -1, 23, 8,  -1, 24, 7,  -1, 25, 6,  -1, 26, 5,  -1, 27, 4,  -1, 28, 3,  -1,
                            29, 2,  -1, 30, 1,  -1, 31, 0,  -1, 48, -1, -1, 49, -1, -1, 50, -1, -1, 51, -1,
                            -1, 52, -1, -1, 53, -1, -1, 54, -1, -1, 55, -1, -1, 56, -1, -1, 57, -1, -1, 58,
                            -1, -1, 59, -1, -1, 60, -1, -1, 61, -1, -1, 62, -1, -1, 63, -1, 47, 32, -1, 46,
                            33, -1, 45, 34, -1, 44, 35, -1, 43, 36, -1, 42, 37, -1, 41, 38, -1, 40, 39, -1 } } },
      /* PS */
      { { 5, 0.5 } }
    };
  } else {
    return new Segmentation{
      5,
      false,
      /* PG */
      { { 1025, 1, 0, 45.7142868, -20 },
        { 1026, 2, 0, 34.2857132, -20 },
        { 1027, 0, 0, 20, -20 },
        { 1033, 1, 0, 5.714285851, -20 },
        { 1034, 2, 0, -5.714285851, -20 },
        { 1035, 0, 0, -20, -20 },
        { 1041, 1, 0, -34.2857132, -20 },
        { 1042, 2, 0, -45.7142868, -20 },
        { 1043, 0, 0, -60, -20 },
        { 1130, 5, 0, 20, -5 },
        { 1131, 3, 0, 28.5714283, -5 },
        { 1132, 4, 0, 40, -5 },
        { 1133, 6, 0, 50, -5 },
        { 1139, 5, 0, -20, -5 },
        { 1140, 3, 0, -11.4285717, -5 },
        { 1141, 4, 0, -8.000008656e-09, -5 },
        { 1142, 6, 0, 10, -5 },
        { 1148, 5, 0, -60, -5 },
        { 1149, 3, 0, -51.42856979, -5 },
        { 1150, 4, 0, -40, -5 },
        { 1151, 6, 0, -30, -5 } },
      /* PGT */
      { /* L1 */ { 20, 4, { 3, 7, 11, 15, 18, 21, 24, 27, 30, 49, 52, 55, 58, 61, 32, 35, 38, 41, 44, 47,
                            2, 6, 10, 14, 17, 20, 23, 26, 29, 48, 51, 54, 57, 60, 63, 34, 37, 40, 43, 46,
                            1, 5, 9,  13, 16, 19, 22, 25, 28, 31, 50, 53, 56, 59, 62, 33, 36, 39, 42, 45,
                            0, 4, 8,  12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 } },
        /* L2 */ { 20, 4, { 2,  5,  8,  11, 14, 17, 20, 23, 26, 29, 48, 51, 54, 57, 60, 63, 35, 39, 43, 47,
                            1,  4,  7,  10, 13, 16, 19, 22, 25, 28, 31, 50, 53, 56, 59, 62, 34, 38, 42, 46,
                            0,  3,  6,  9,  12, 15, 18, 21, 24, 27, 30, 49, 52, 55, 58, 61, 33, 37, 41, 45,
                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 32, 36, 40, 44 } },
        /* O3 */ { 16, 4, { 3,  7,  11, 15, 19, 23, 27, 31, 51, 55, 59, 63, 35, 39, 43, 47, 2,  6,  10, 14, 18, 22,
                            26, 30, 50, 54, 58, 62, 34, 38, 42, 46, 1,  5,  9,  13, 17, 21, 25, 29, 49, 53, 57, 61,
                            33, 37, 41, 45, 0,  4,  8,  12, 16, 20, 24, 28, 48, 52, 56, 60, 32, 36, 40, 44 } },
        /* P1 */ { 16, 5, { 47, 46, 41, 36, 63, 58, 53, 48, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 42, 37,
                            32, 59, 54, 49, 28, 24, 20, 16, 12, 8,  4,  0,  -1, -1, 43, 38, 33, 60, 55, 50,
                            29, 25, 21, 17, 13, 9,  5,  1,  -1, -1, 44, 39, 34, 61, 56, 51, 30, 26, 22, 18,
                            14, 10, 6,  2,  -1, -1, 45, 40, 35, 62, 57, 52, 31, 27, 23, 19, 15, 11, 7,  3 } },
        /* P2 */ { 16, 5, { -1, -1, -1, -1, -1, -1, -1, -1, 27, 22, 17, 12, 7,  2,  1,  0,  44, 40, 36, 32,
                            60, 56, 52, 48, 28, 23, 18, 13, 8,  3,  -1, -1, 45, 41, 37, 33, 61, 57, 53, 49,
                            29, 24, 19, 14, 9,  4,  -1, -1, 46, 42, 38, 34, 62, 58, 54, 50, 30, 25, 20, 15,
                            10, 5,  -1, -1, 47, 43, 39, 35, 63, 59, 55, 51, 31, 26, 21, 16, 11, 6,  -1, -1 } },
        /* Q1 */ { 14, 5, { -1, -1, -1, -1, 59, 54, 49, 28, 23, 18, 13, 8,  -1, -1, 44, 40, 36, 32,
                            60, 55, 50, 29, 24, 19, 14, 9,  4,  0,  45, 41, 37, 33, 61, 56, 51, 30,
                            25, 20, 15, 10, 5,  1,  46, 42, 38, 34, 62, 57, 52, 31, 26, 21, 16, 11,
                            6,  2,  47, 43, 39, 35, 63, 58, 53, 48, 27, 22, 17, 12, 7,  3 } },
        /* Q2 */ { 14, 5, { -1, -1, 35, 62, 57, 52, 31, 26, 21, 16, -1, -1, -1, -1, 44, 40, 36, 63,
                            58, 53, 48, 27, 22, 17, 12, 8,  4,  0,  45, 41, 37, 32, 59, 54, 49, 28,
                            23, 18, 13, 9,  5,  1,  46, 42, 38, 33, 60, 55, 50, 29, 24, 19, 14, 10,
                            6,  2,  47, 43, 39, 34, 61, 56, 51, 30, 25, 20, 15, 11, 7,  3 } } },
      /* PS */
      { { 0.7142857313, 5 } }
    };
  }
}
class SegmentationCreatorRegisterCreateSegType5
{
 public:
  SegmentationCreatorRegisterCreateSegType5() { registerSegmentationCreator(5, createSegType5); }
} aSegmentationCreatorRegisterCreateSegType5;

} // namespace impl2
} // namespace mapping
} // namespace mch
} // namespace o2
