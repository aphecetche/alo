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
Segmentation* createSegType0(bool isBendingPlane)
{
  if (isBendingPlane) {
    return new Segmentation{
      0,
      true,
      /* PG */
      { { 1, 6, 3, 80.64, 0 },       { 1, 9, 1, 88.2, 0 },        { 2, 3, 3, 70.56, 0 },
        { 3, 3, 2, 65.52, 0 },       { 4, 3, 2, 60.48, 0 },       { 5, 3, 2, 55.44, 0 },
        { 6, 3, 0, 52.92, 0 },       { 7, 3, 0, 50.4, 0 },        { 8, 3, 0, 47.88, 0 },
        { 9, 3, 0, 45.36, 0 },       { 10, 3, 0, 42.84, 0 },      { 11, 3, 0, 40.32, 0 },
        { 12, 3, 0, 37.8, 0 },       { 13, 3, 0, 35.28, 0 },      { 14, 3, 0, 32.76, 0 },
        { 15, 3, 0, 30.24, 0 },      { 16, 3, 0, 27.72, 0 },      { 17, 3, 0, 25.2, 0 },
        { 18, 3, 0, 22.68, 0 },      { 19, 3, 0, 20.16, 0 },      { 20, 3, 0, 17.64, 0 },
        { 27, 6, 3, 80.64, 6.72 },   { 27, 9, 1, 88.2, 6.72 },    { 28, 3, 3, 70.56, 6.72 },
        { 29, 3, 2, 65.52, 6.72 },   { 30, 3, 2, 60.48, 6.72 },   { 31, 3, 2, 55.44, 6.72 },
        { 32, 3, 0, 52.92, 6.72 },   { 33, 3, 0, 50.4, 6.72 },    { 34, 3, 0, 47.88, 6.72 },
        { 35, 3, 0, 45.36, 6.72 },   { 36, 3, 0, 42.84, 6.72 },   { 37, 3, 0, 40.32, 6.72 },
        { 38, 3, 0, 37.8, 6.72 },    { 39, 3, 0, 35.28, 6.72 },   { 40, 3, 0, 32.76, 6.72 },
        { 41, 3, 0, 30.24, 6.72 },   { 42, 3, 0, 27.72, 6.72 },   { 43, 3, 0, 25.2, 6.72 },
        { 44, 3, 0, 22.68, 6.72 },   { 45, 3, 0, 20.16, 6.72 },   { 46, 3, 0, 17.64, 6.72 },
        { 47, 2, 0, 13.86, 6.72 },   { 53, 6, 3, 80.64, 13.44 },  { 53, 9, 1, 88.2, 13.44 },
        { 54, 3, 3, 70.56, 13.44 },  { 55, 3, 2, 65.52, 13.44 },  { 56, 3, 2, 60.48, 13.44 },
        { 57, 3, 2, 55.44, 13.44 },  { 58, 3, 0, 52.92, 13.44 },  { 59, 3, 0, 50.4, 13.44 },
        { 60, 3, 0, 47.88, 13.44 },  { 61, 3, 0, 45.36, 13.44 },  { 62, 3, 0, 42.84, 13.44 },
        { 63, 3, 0, 40.32, 13.44 },  { 64, 3, 0, 37.8, 13.44 },   { 65, 3, 0, 35.28, 13.44 },
        { 66, 3, 0, 32.76, 13.44 },  { 67, 3, 0, 30.24, 13.44 },  { 68, 3, 0, 27.72, 13.44 },
        { 69, 3, 0, 25.2, 13.44 },   { 70, 3, 0, 22.68, 13.44 },  { 71, 3, 0, 20.16, 13.44 },
        { 72, 3, 0, 17.64, 13.44 },  { 73, 3, 0, 15.12, 13.44 },  { 74, 3, 0, 12.6, 13.44 },
        { 75, 1, 0, 9.45, 13.44 },   { 76, 0, 0, 0, 17.22 },      { 79, 6, 3, 80.64, 20.16 },
        { 79, 9, 1, 88.2, 20.16 },   { 80, 3, 3, 70.56, 20.16 },  { 81, 3, 2, 65.52, 20.16 },
        { 82, 3, 2, 60.48, 20.16 },  { 83, 3, 2, 55.44, 20.16 },  { 84, 3, 2, 50.4, 20.16 },
        { 85, 3, 0, 47.88, 20.16 },  { 86, 3, 0, 45.36, 20.16 },  { 87, 3, 0, 42.84, 20.16 },
        { 88, 3, 0, 40.32, 20.16 },  { 89, 3, 0, 37.8, 20.16 },   { 90, 3, 0, 35.28, 20.16 },
        { 91, 3, 0, 32.76, 20.16 },  { 92, 3, 0, 30.24, 20.16 },  { 93, 3, 0, 27.72, 20.16 },
        { 94, 3, 0, 25.2, 20.16 },   { 95, 3, 0, 22.68, 20.16 },  { 96, 3, 0, 20.16, 20.16 },
        { 97, 3, 0, 17.64, 20.16 },  { 98, 3, 0, 15.12, 20.16 },  { 99, 3, 0, 12.6, 20.16 },
        { 100, 3, 0, 10.08, 20.16 }, { 101, 3, 0, 7.56, 20.16 },  { 102, 3, 0, 5.04, 20.16 },
        { 103, 3, 0, 2.52, 20.16 },  { 104, 3, 0, 0, 20.16 },     { 105, 6, 3, 80.64, 26.88 },
        { 105, 9, 1, 88.2, 26.88 },  { 106, 3, 3, 70.56, 26.88 }, { 107, 3, 2, 65.52, 26.88 },
        { 108, 3, 2, 60.48, 26.88 }, { 109, 3, 2, 55.44, 26.88 }, { 110, 3, 2, 50.4, 26.88 },
        { 111, 3, 0, 47.88, 26.88 }, { 112, 3, 0, 45.36, 26.88 }, { 113, 3, 0, 42.84, 26.88 },
        { 114, 3, 0, 40.32, 26.88 }, { 115, 3, 0, 37.8, 26.88 },  { 116, 3, 0, 35.28, 26.88 },
        { 117, 3, 0, 32.76, 26.88 }, { 118, 3, 0, 30.24, 26.88 }, { 119, 3, 0, 27.72, 26.88 },
        { 120, 3, 0, 25.2, 26.88 },  { 121, 3, 0, 22.68, 26.88 }, { 122, 3, 0, 20.16, 26.88 },
        { 123, 3, 0, 17.64, 26.88 }, { 124, 3, 0, 15.12, 26.88 }, { 125, 3, 0, 12.6, 26.88 },
        { 126, 3, 0, 10.08, 26.88 }, { 127, 3, 0, 7.56, 26.88 },  { 128, 3, 0, 5.04, 26.88 },
        { 129, 3, 0, 2.52, 26.88 },  { 130, 3, 0, 0, 26.88 },     { 131, 6, 3, 80.64, 33.6 },
        { 131, 9, 1, 88.2, 33.6 },   { 132, 3, 3, 70.56, 33.6 },  { 133, 3, 3, 60.48, 33.6 },
        { 134, 3, 2, 55.44, 33.6 },  { 135, 3, 2, 50.4, 33.6 },   { 136, 3, 2, 45.36, 33.6 },
        { 137, 3, 0, 42.84, 33.6 },  { 138, 3, 0, 40.32, 33.6 },  { 139, 3, 0, 37.8, 33.6 },
        { 140, 3, 0, 35.28, 33.6 },  { 141, 3, 0, 32.76, 33.6 },  { 142, 3, 0, 30.24, 33.6 },
        { 143, 3, 0, 27.72, 33.6 },  { 144, 3, 0, 25.2, 33.6 },   { 145, 3, 0, 22.68, 33.6 },
        { 146, 3, 0, 20.16, 33.6 },  { 147, 3, 0, 17.64, 33.6 },  { 148, 3, 0, 15.12, 33.6 },
        { 149, 3, 0, 12.6, 33.6 },   { 150, 3, 0, 10.08, 33.6 },  { 151, 3, 0, 7.56, 33.6 },
        { 152, 3, 0, 5.04, 33.6 },   { 153, 3, 0, 2.52, 33.6 },   { 154, 3, 0, 0, 33.6 },
        { 157, 6, 3, 80.64, 40.32 }, { 157, 9, 1, 88.2, 40.32 },  { 158, 3, 3, 70.56, 40.32 },
        { 159, 3, 3, 60.48, 40.32 }, { 160, 3, 2, 55.44, 40.32 }, { 161, 3, 2, 50.4, 40.32 },
        { 162, 3, 2, 45.36, 40.32 }, { 163, 3, 0, 42.84, 40.32 }, { 164, 3, 0, 40.32, 40.32 },
        { 165, 3, 0, 37.8, 40.32 },  { 166, 3, 0, 35.28, 40.32 }, { 167, 3, 0, 32.76, 40.32 },
        { 168, 3, 0, 30.24, 40.32 }, { 169, 3, 0, 27.72, 40.32 }, { 170, 3, 0, 25.2, 40.32 },
        { 171, 3, 0, 22.68, 40.32 }, { 172, 3, 0, 20.16, 40.32 }, { 173, 3, 0, 17.64, 40.32 },
        { 174, 3, 0, 15.12, 40.32 }, { 175, 3, 0, 12.6, 40.32 },  { 176, 3, 0, 10.08, 40.32 },
        { 177, 3, 0, 7.56, 40.32 },  { 178, 3, 0, 5.04, 40.32 },  { 179, 3, 0, 2.52, 40.32 },
        { 180, 3, 0, 0, 40.32 },     { 183, 3, 3, 70.56, 47.04 }, { 184, 3, 3, 60.48, 47.04 },
        { 185, 3, 3, 50.4, 47.04 },  { 186, 3, 2, 45.36, 47.04 }, { 187, 3, 2, 40.32, 47.04 },
        { 188, 3, 2, 35.28, 47.04 }, { 189, 3, 2, 30.24, 47.04 }, { 190, 3, 2, 25.2, 47.04 },
        { 191, 3, 0, 22.68, 47.04 }, { 192, 3, 0, 20.16, 47.04 }, { 193, 3, 0, 17.64, 47.04 },
        { 194, 3, 0, 15.12, 47.04 }, { 195, 3, 0, 12.6, 47.04 },  { 196, 3, 0, 10.08, 47.04 },
        { 197, 3, 0, 7.56, 47.04 },  { 198, 3, 0, 5.04, 47.04 },  { 199, 3, 0, 2.52, 47.04 },
        { 200, 3, 0, 0, 47.04 },     { 201, 8, 3, 70.56, 53.76 }, { 202, 3, 3, 60.48, 53.76 },
        { 203, 3, 3, 50.4, 53.76 },  { 204, 3, 2, 45.36, 53.76 }, { 205, 3, 2, 40.32, 53.76 },
        { 206, 3, 2, 35.28, 53.76 }, { 207, 3, 2, 30.24, 53.76 }, { 208, 3, 2, 25.2, 53.76 },
        { 209, 3, 2, 20.16, 53.76 }, { 210, 3, 2, 15.12, 53.76 }, { 211, 3, 2, 10.08, 53.76 },
        { 212, 3, 2, 5.04, 53.76 },  { 213, 3, 2, 0, 53.76 },     { 214, 4, 3, 60.48, 60.48 },
        { 215, 3, 3, 50.4, 60.48 },  { 216, 3, 3, 40.32, 60.48 }, { 217, 3, 3, 30.24, 60.48 },
        { 218, 3, 2, 25.2, 60.48 },  { 219, 3, 2, 20.16, 60.48 }, { 220, 3, 2, 15.12, 60.48 },
        { 221, 3, 2, 10.08, 60.48 }, { 222, 3, 2, 5.04, 60.48 },  { 223, 3, 2, 0, 60.48 },
        { 224, 3, 3, 50.4, 67.2 },   { 225, 3, 3, 40.32, 67.2 },  { 226, 3, 3, 30.24, 67.2 },
        { 227, 3, 3, 20.16, 67.2 },  { 228, 3, 2, 15.12, 67.2 },  { 229, 3, 2, 10.08, 67.2 },
        { 230, 3, 2, 5.04, 67.2 },   { 231, 3, 2, 0, 67.2 },      { 232, 7, 3, 50.4, 73.92 },
        { 233, 3, 3, 40.32, 73.92 }, { 234, 3, 3, 30.24, 73.92 }, { 235, 3, 3, 20.16, 73.92 },
        { 236, 3, 3, 10.08, 73.92 }, { 237, 3, 3, 0, 73.92 },     { 238, 5, 3, 37.8, 80.64 },
        { 239, 5, 3, 30.24, 80.64 }, { 240, 5, 3, 22.68, 80.64 }, { 241, 5, 3, 15.12, 80.64 },
        { 242, 5, 3, 7.56, 80.64 },  { 243, 5, 3, 0, 80.64 } },
      /* PGT */
      { /* 1BA */ { 12, 7, { -1, -1, -1, -1, -1, -1, -1, -1, -1, 26, 24, -1, -1, -1, -1, -1, -1, -1, -1, 30, 27,
                             52, 57, 18, -1, -1, -1, 7,  4,  1,  46, 48, 53, 56, 22, 62, 39, 38, 9,  8,  41, 45,
                             0,  50, 54, 21, 60, 63, 10, 35, 13, 36, 6,  44, 47, 29, 25, 59, 17, 16, 12, 32, 34,
                             11, 42, 3,  31, 51, 23, 20, 61, 19, 15, 33, 14, 37, 40, 5,  43, 2,  49, 28, 55, 58 } },
        /* 1BB */ { 5, 16, { -1, -1, -1, -1, 63, -1, -1, -1, 17, 16, -1, -1, -1, 62, 61, -1, -1, 18, 19, 58,
                             -1, 60, 59, 20, 22, -1, 21, 57, 56, 55, -1, 23, 24, 54, 25, -1, 53, 26, 27, 52,
                             51, 50, 29, 30, 28, 49, 48, 31, 0,  2,  47, 46, 1,  44, 43, 45, 3,  4,  41, 5,
                             42, 6,  40, 39, 8,  7,  38, 9,  37, 11, 10, 36, 35, 12, 34, 13, 33, 32, 15, 14 } },
        /* 1BC */ { 6, 16, { -1, -1, -1, -1, 63, 16, -1, -1, -1, -1, 17, 61, -1, -1, -1, -1, 62, 19, -1, -1,
                             -1, 18, 60, 58, -1, -1, -1, 59, 20, 22, -1, -1, -1, 21, 57, 55, -1, -1, 56, 23,
                             24, 25, -1, -1, 54, 53, 26, 52, -1, -1, 27, 51, 29, 28, -1, -1, 50, 49, 47, 31,
                             -1, 30, 48, 46, 44, 2,  -1, 0,  1,  3,  4,  5,  -1, 45, 43, 41, 39, 40, 42, 6,
                             7,  38, 36, 8,  10, 11, 35, 12, 34, 37, 13, 33, 14, 32, 15, 9 } },
        /* 1BD */ { 4, 16, { 18, 63, 16, 20, 57, 62, 61, 58, 21, 17, 60, 23, 24, 59, 19, 55, 53, 56, 22, 26, 50, 54,
                             25, 51, 30, 27, 52, 28, 47, 48, 29, 49, 1,  0,  31, 46, 44, 3,  45, 2,  41, 4,  42, 43,
                             39, 7,  40, 5,  10, 38, 8,  6,  36, 35, 9,  37, 12, 13, 34, 11, 33, 32, 15, 14 } },
        /* 1BE */ { 3, 16, { 18, 63, 16, 57, 62, 61, 21, 17, 60, 24, 59, 19, 53, 56, 22, 50,
                             54, 25, 30, 27, 52, 47, 48, 29, 1,  0,  31, 44, 3,  45, 41, 4,
                             42, 39, 7,  40, 10, 38, 8,  36, 35, 9,  12, 13, 34, 33, 32, 15 } },
        /* 1BF */ { 3, 21, { 59, 63, 19, 57, 17, 20, 21, 16, 22, 54, 62, 23, 24, 61, 55, 53, 18, 26, 27, 58, 52,
                             50, 56, 28, 30, 25, 49, 48, 51, 31, 0,  29, 46, 1,  47, 45, 44, 2,  43, 41, 4,  3,
                             6,  42, 5,  7,  38, 40, 9,  37, 8,  10, 35, 39, 36, 14, 11, 13, 15, 34, 32, 33, 12 } },
        /* 1BG */ { 3, 16, { 18, 63, 16, 57, 62, 61, 21, 17, 60, 24, 59, 19, 53, 56, 22, 50,
                             54, 25, 30, 27, 52, 47, 48, 29, 1,  0,  31, 44, 3,  45, 41, 4,
                             42, 39, 7,  40, 10, 38, 8,  36, 35, 9,  12, 13, 34, 33, 32, 15 } },
        /* 1BH */ { 3, 16, { 62, 63, 16, 18, 17, 19, 59, 60, 58, 21, 22, 23, 56, 54, 25, 24,
                             26, 52, 27, 28, 29, 50, 48, 31, 30, 46, 1,  47, 2,  3,  44, 4,
                             42, 41, 40, 7,  38, 9,  8,  10, 11, 36, 32, 34, -1, 15, 33, -1 } },
        /* 1BI */ { 2, 16, { 62, 63, 18, 17, 59, 60, 21, 22, 56, 54, 24, 26, 27, 28, 50, 48,
                             30, 46, 47, 2,  44, 4,  41, 40, 38, 9,  10, 11, 32, 34, 15, 33 } },
        /* 1BG */ { 1, 16, { 20, 58, 23, 55, 26, 51, 28, 49, 46, 2, 43, 5, 6, 37, 11, 14 } } },
      /* PS */
      { { 0.63, 0.42 }, { 0.84, 0.42 }, { 1.26, 0.42 }, { 2.52, 0.42 } }
    };
  } else {
    return new Segmentation{
      0,
      false,
      /* PG */
      { { 1025, 8, 2, 80.325, 0.21 },   { 1026, 5, 2, 70.245, 0.21 },   { 1027, 4, 1, 65.205, 0.21 },
        { 1028, 4, 1, 60.165, 0.21 },   { 1029, 4, 1, 55.125, 0.21 },   { 1030, 0, 0, 52.605, 0.21 },
        { 1031, 0, 0, 50.085, 0.21 },   { 1032, 0, 0, 47.565, 0.21 },   { 1033, 0, 0, 45.045, 0.21 },
        { 1034, 0, 0, 42.525, 0.21 },   { 1035, 0, 0, 40.005, 0.21 },   { 1036, 0, 0, 37.485, 0.21 },
        { 1037, 0, 0, 34.965, 0.21 },   { 1038, 0, 0, 32.445, 0.21 },   { 1039, 0, 0, 29.925, 0.21 },
        { 1040, 0, 0, 27.405, 0.21 },   { 1041, 0, 0, 24.885, 0.21 },   { 1042, 0, 0, 22.365, 0.21 },
        { 1043, 0, 0, 19.845, 0.21 },   { 1044, 11, 0, 16.695, 0.21 },  { 1051, 8, 2, 80.325, 6.93 },
        { 1052, 5, 2, 70.245, 6.93 },   { 1053, 4, 1, 65.205, 6.93 },   { 1054, 4, 1, 60.165, 6.93 },
        { 1055, 4, 1, 55.125, 6.93 },   { 1056, 0, 0, 52.605, 6.93 },   { 1057, 0, 0, 50.085, 6.93 },
        { 1058, 0, 0, 47.565, 6.93 },   { 1059, 0, 0, 45.045, 6.93 },   { 1060, 0, 0, 42.525, 6.93 },
        { 1061, 0, 0, 40.005, 6.93 },   { 1062, 0, 0, 37.485, 6.93 },   { 1063, 0, 0, 34.965, 6.93 },
        { 1064, 0, 0, 32.445, 6.93 },   { 1065, 0, 0, 29.925, 6.93 },   { 1066, 0, 0, 27.405, 6.93 },
        { 1067, 0, 0, 24.885, 6.93 },   { 1068, 0, 0, 22.365, 6.93 },   { 1069, 0, 0, 19.845, 6.93 },
        { 1070, 12, 0, 16.695, 7.35 },  { 1071, 3, 0, 11.655, 7.35 },   { 1077, 8, 2, 80.325, 13.65 },
        { 1078, 5, 2, 70.245, 13.65 },  { 1079, 4, 1, 65.205, 13.65 },  { 1080, 4, 1, 60.165, 13.65 },
        { 1081, 4, 1, 55.125, 13.65 },  { 1082, 0, 0, 52.605, 13.65 },  { 1083, 0, 0, 50.085, 13.65 },
        { 1084, 0, 0, 47.565, 13.65 },  { 1085, 0, 0, 45.045, 13.65 },  { 1086, 0, 0, 42.525, 13.65 },
        { 1087, 0, 0, 40.005, 13.65 },  { 1088, 0, 0, 37.485, 13.65 },  { 1089, 0, 0, 34.965, 13.65 },
        { 1090, 0, 0, 32.445, 13.65 },  { 1091, 0, 0, 29.925, 13.65 },  { 1092, 0, 0, 27.405, 13.65 },
        { 1093, 0, 0, 24.885, 13.65 },  { 1094, 0, 0, 22.365, 13.65 },  { 1095, 0, 0, 19.845, 13.65 },
        { 1096, 0, 0, 17.325, 13.65 },  { 1097, 0, 0, 14.805, 13.65 },  { 1098, 13, 0, 11.655, 14.49 },
        { 1099, 2, 0, 5.985, 14.49 },   { 1100, 1, 0, 1.575, 17.85 },   { 1103, 8, 2, 80.325, 20.37 },
        { 1104, 5, 2, 70.245, 20.37 },  { 1105, 4, 1, 65.205, 20.37 },  { 1106, 4, 1, 60.165, 20.37 },
        { 1107, 4, 1, 55.125, 20.37 },  { 1108, 4, 1, 50.085, 20.37 },  { 1109, 0, 0, 47.565, 20.37 },
        { 1110, 0, 0, 45.045, 20.37 },  { 1111, 0, 0, 42.525, 20.37 },  { 1112, 0, 0, 40.005, 20.37 },
        { 1113, 0, 0, 37.485, 20.37 },  { 1114, 0, 0, 34.965, 20.37 },  { 1115, 0, 0, 32.445, 20.37 },
        { 1116, 0, 0, 29.925, 20.37 },  { 1117, 0, 0, 27.405, 20.37 },  { 1118, 0, 0, 24.885, 20.37 },
        { 1119, 0, 0, 22.365, 20.37 },  { 1120, 0, 0, 19.845, 20.37 },  { 1121, 0, 0, 17.325, 20.37 },
        { 1122, 0, 0, 14.805, 20.37 },  { 1123, 0, 0, 12.285, 20.37 },  { 1124, 0, 0, 9.765, 20.37 },
        { 1125, 0, 0, 7.245, 20.37 },   { 1126, 0, 0, 4.725, 20.37 },   { 1127, 0, 0, 2.205, 20.37 },
        { 1128, 0, 0, -0.315, 20.37 },  { 1129, 8, 2, 80.325, 27.09 },  { 1130, 5, 2, 70.245, 27.09 },
        { 1131, 4, 1, 65.205, 27.09 },  { 1132, 4, 1, 60.165, 27.09 },  { 1133, 4, 1, 55.125, 27.09 },
        { 1134, 4, 1, 50.085, 27.09 },  { 1135, 4, 1, 45.045, 27.09 },  { 1136, 0, 0, 42.525, 27.09 },
        { 1137, 0, 0, 40.005, 27.09 },  { 1138, 0, 0, 37.485, 27.09 },  { 1139, 0, 0, 34.965, 27.09 },
        { 1140, 0, 0, 32.445, 27.09 },  { 1141, 0, 0, 29.925, 27.09 },  { 1142, 0, 0, 27.405, 27.09 },
        { 1143, 0, 0, 24.885, 27.09 },  { 1144, 0, 0, 22.365, 27.09 },  { 1145, 0, 0, 19.845, 27.09 },
        { 1146, 0, 0, 17.325, 27.09 },  { 1147, 0, 0, 14.805, 27.09 },  { 1148, 0, 0, 12.285, 27.09 },
        { 1149, 0, 0, 9.765, 27.09 },   { 1150, 0, 0, 7.245, 27.09 },   { 1151, 0, 0, 4.725, 27.09 },
        { 1152, 0, 0, 2.205, 27.09 },   { 1153, 0, 0, -0.315, 27.09 },  { 1155, 8, 2, 80.325, 33.81 },
        { 1156, 5, 2, 70.245, 33.81 },  { 1157, 5, 2, 60.165, 33.81 },  { 1158, 4, 1, 55.125, 33.81 },
        { 1159, 4, 1, 50.085, 33.81 },  { 1160, 4, 1, 45.045, 33.81 },  { 1161, 0, 0, 42.525, 33.81 },
        { 1162, 0, 0, 40.005, 33.81 },  { 1163, 0, 0, 37.485, 33.81 },  { 1164, 0, 0, 34.965, 33.81 },
        { 1165, 0, 0, 32.445, 33.81 },  { 1166, 0, 0, 29.925, 33.81 },  { 1167, 0, 0, 27.405, 33.81 },
        { 1168, 0, 0, 24.885, 33.81 },  { 1169, 0, 0, 22.365, 33.81 },  { 1170, 0, 0, 19.845, 33.81 },
        { 1171, 0, 0, 17.325, 33.81 },  { 1172, 0, 0, 14.805, 33.81 },  { 1173, 0, 0, 12.285, 33.81 },
        { 1174, 0, 0, 9.765, 33.81 },   { 1175, 0, 0, 7.245, 33.81 },   { 1176, 0, 0, 4.725, 33.81 },
        { 1177, 0, 0, 2.205, 33.81 },   { 1178, 0, 0, -0.315, 33.81 },  { 1181, 8, 2, 80.325, 40.53 },
        { 1182, 5, 2, 70.245, 40.53 },  { 1183, 5, 2, 60.165, 40.53 },  { 1184, 4, 1, 55.125, 40.53 },
        { 1185, 4, 1, 50.085, 40.53 },  { 1186, 4, 1, 45.045, 40.53 },  { 1187, 0, 0, 42.525, 40.53 },
        { 1188, 0, 0, 40.005, 40.53 },  { 1189, 0, 0, 37.485, 40.53 },  { 1190, 0, 0, 34.965, 40.53 },
        { 1191, 0, 0, 32.445, 40.53 },  { 1192, 0, 0, 29.925, 40.53 },  { 1193, 0, 0, 27.405, 40.53 },
        { 1194, 0, 0, 24.885, 40.53 },  { 1195, 0, 0, 22.365, 40.53 },  { 1196, 0, 0, 19.845, 40.53 },
        { 1197, 0, 0, 17.325, 40.53 },  { 1198, 0, 0, 14.805, 40.53 },  { 1199, 0, 0, 12.285, 40.53 },
        { 1200, 0, 0, 9.765, 40.53 },   { 1201, 0, 0, 7.245, 40.53 },   { 1202, 0, 0, 4.725, 40.53 },
        { 1203, 0, 0, 2.205, 40.53 },   { 1204, 0, 0, -0.315, 40.53 },  { 1207, 5, 2, 70.245, 47.25 },
        { 1208, 5, 2, 60.165, 47.25 },  { 1209, 5, 2, 50.085, 47.25 },  { 1210, 4, 1, 45.045, 47.25 },
        { 1211, 4, 1, 40.005, 47.25 },  { 1212, 4, 1, 34.965, 47.25 },  { 1213, 4, 1, 29.925, 47.25 },
        { 1214, 4, 1, 24.885, 47.25 },  { 1215, 0, 0, 22.365, 47.25 },  { 1216, 0, 0, 19.845, 47.25 },
        { 1217, 0, 0, 17.325, 47.25 },  { 1218, 0, 0, 14.805, 47.25 },  { 1219, 0, 0, 12.285, 47.25 },
        { 1220, 0, 0, 9.765, 47.25 },   { 1221, 0, 0, 7.245, 47.25 },   { 1222, 0, 0, 4.725, 47.25 },
        { 1223, 0, 0, 2.205, 47.25 },   { 1224, 0, 0, -0.315, 47.25 },  { 1225, 10, 2, 70.245, 53.97 },
        { 1226, 5, 2, 60.165, 53.97 },  { 1227, 5, 2, 50.085, 53.97 },  { 1228, 4, 1, 45.045, 53.97 },
        { 1229, 4, 1, 40.005, 53.97 },  { 1230, 4, 1, 34.965, 53.97 },  { 1231, 4, 1, 29.925, 53.97 },
        { 1232, 4, 1, 24.885, 53.97 },  { 1233, 4, 1, 19.845, 53.97 },  { 1234, 4, 1, 14.805, 53.97 },
        { 1235, 4, 1, 9.765, 53.97 },   { 1236, 4, 1, 4.725, 53.97 },   { 1237, 4, 1, -0.315, 53.97 },
        { 1238, 5, 2, 60.165, 60.69 },  { 1239, 5, 2, 50.085, 60.69 },  { 1240, 5, 2, 40.005, 60.69 },
        { 1241, 4, 1, 34.965, 60.69 },  { 1242, 4, 1, 29.925, 60.69 },  { 1243, 4, 1, 24.885, 60.69 },
        { 1244, 4, 1, 19.845, 60.69 },  { 1245, 4, 1, 14.805, 60.69 },  { 1246, 4, 1, 9.765, 60.69 },
        { 1247, 4, 1, 4.725, 60.69 },   { 1248, 4, 1, -0.315, 60.69 },  { 1249, 5, 2, 50.085, 67.41 },
        { 1250, 5, 2, 40.005, 67.41 },  { 1251, 5, 2, 29.925, 67.41 },  { 1252, 5, 2, 19.845, 67.41 },
        { 1253, 4, 1, 14.805, 67.41 },  { 1254, 4, 1, 9.765, 67.41 },   { 1255, 4, 1, 4.725, 67.41 },
        { 1256, 4, 1, -0.315, 67.41 },  { 1257, 9, 2, 50.085, 74.13 },  { 1258, 5, 2, 40.005, 74.13 },
        { 1259, 5, 2, 29.925, 74.13 },  { 1260, 5, 2, 19.845, 74.13 },  { 1261, 5, 2, 9.765, 74.13 },
        { 1262, 5, 2, -0.315, 74.13 },  { 1263, 7, 2, 40.005, 80.85 },  { 1263, 15, 3, 40.005, 85.89 },
        { 1264, 6, 2, 29.925, 80.85 },  { 1264, 14, 3, 29.925, 85.89 }, { 1265, 6, 2, 19.845, 80.85 },
        { 1265, 14, 3, 19.845, 85.89 }, { 1266, 6, 2, 9.765, 80.85 },   { 1266, 14, 3, 9.765, 85.89 },
        { 1267, 6, 2, -0.315, 80.85 },  { 1267, 14, 3, -0.315, 85.89 } },
      /* PGT */
      { /* 1NA */ { 4, 16, { 20, 16, 63, 18, 58, 61, 62, 57, 23, 60, 17, 21, 55, 19, 59, 24, 52, 22, 56, 53, 26, 25,
                             54, 27, 49, 51, 50, 47, 28, 48, 29, 30, 31, 46, 45, 44, 43, 0,  42, 1,  2,  3,  41, 4,
                             40, 5,  38, 6,  39, 37, 7,  35, 8,  36, 10, 9,  34, 11, 32, 33, 13, 15, 14, 12 } },
        /* 1NB */ { 13, 6, { -1, -1, -1, -1, -1, 51, 54, -1, -1, -1, -1, -1, -1, -1, -1, 45, 48, 27, 31, 28,
                             49, 55, 22, 19, 61, 17, -1, 43, 2,  42, 40, 3,  30, 46, 26, 25, 58, 60, 16, 37,
                             6,  15, 13, 12, 8,  5,  0,  29, 52, 23, 20, 63, 34, 11, 33, 14, 32, 9,  39, 1,
                             47, 50, 56, 59, 62, 10, 36, 35, 38, 7,  41, 4,  44, 53, 24, 57, 21, 18 } },
        /* 1NC */ { 10, 14, { -1, -1, -1, -1, -1, -1, -1, -1, 63, -1, -1, -1, -1, -1, -1, -1, -1, 16, 62, -1,
                              -1, -1, -1, -1, -1, -1, 61, 17, 19, -1, -1, -1, -1, -1, -1, 60, 18, 20, 59, -1,
                              -1, -1, -1, -1, 58, 22, 57, 21, 56, -1, -1, -1, -1, 23, 25, 54, 24, 51, 52, -1,
                              -1, 55, 53, 49, 26, 50, 27, 48, 28, -1, 47, 29, 46, 30, 0,  43, 44, 31, 45, -1,
                              42, 41, 40, 39, 5,  6,  3,  2,  1,  -1, -1, -1, -1, -1, -1, -1, 37, 38, 4,  -1,
                              -1, -1, -1, -1, -1, -1, 8,  36, 7,  -1, -1, -1, -1, -1, -1, -1, 34, 9,  35, -1,
                              -1, -1, -1, -1, -1, -1, 10, 33, 11, 32, -1, -1, -1, -1, -1, -1, 12, 13, 15, 14 } },
        /* 1ND */ { 9, 21, { -1, -1, -1, -1, -1, -1, -1, -1, 63, -1, -1, -1, -1, -1, -1, -1, -1, 62, -1, -1, -1,
                             -1, -1, -1, -1, 16, 17, -1, -1, -1, -1, -1, -1, -1, 61, 18, -1, -1, -1, -1, -1, -1,
                             -1, 60, 59, -1, -1, -1, -1, -1, -1, 19, 57, 21, -1, -1, -1, -1, -1, -1, 55, 54, 24,
                             -1, -1, -1, -1, -1, 58, 51, 50, 27, -1, -1, -1, -1, -1, 56, 28, 46, -1, -1, -1, -1,
                             -1, 20, 53, 47, 43, -1, -1, -1, -1, -1, 25, 29, 0,  42, -1, -1, -1, -1, 22, 26, 30,
                             40, 39, -1, -1, -1, -1, 52, 45, 3,  37, 7,  -1, -1, -1, 23, 44, 41, 8,  35, 33, -1,
                             -1, -1, 2,  5,  36, 10, 14, 15, -1, -1, 48, 1,  6,  34, -1, -1, -1, -1, 49, 31, 4,
                             38, 32, -1, -1, -1, -1, -1, -1, -1, -1, 12, -1, -1, -1, -1, -1, -1, -1, -1, 13, -1,
                             -1, -1, -1, -1, -1, -1, -1, 11, -1, -1, -1, -1, -1, -1, -1, -1, 9,  -1, -1, -1, -1 } },
        /* 1NE */ { 8, 8, { 61, 17, 16, 63, 62, 60, 59, 57, 22, 58, 20, 19, 18, 21, 23, 56, 49, 51, 25, 55, 54, 24,
                            52, 53, 46, 29, 28, 48, 50, 26, 27, 47, 42, 43, 0,  31, 30, 45, 44, 1,  5,  40, 3,  2,
                            41, 4,  39, 38, 6,  37, 8,  36, 7,  35, 9,  10, 34, 33, 11, 13, 32, 12, 14, 15 } },
        /* 1NF */ { 16, 4, { 22, 57, 58, 19, 61, 17, 16, 63, 62, 60, 18, 59, 20, 21, 56, 24, 28, 49, 26, 51, 52, 25,
                             55, 23, 54, 53, 50, 48, 27, 47, 30, 44, 40, 2,  43, 29, 46, 45, 31, 0,  1,  42, 41, 3,
                             4,  39, 5,  38, 13, 11, 6,  37, 8,  9,  34, 33, 35, 36, 7,  10, 32, 12, 15, 14 } },
        /* 1NG */ { 16, 3, { 54, 23, 22, 19, 61, 17, 16, 63, 62, 60, 18, 20, 59, 21, 56, 24,
                             28, 49, 26, 51, 52, 25, 55, 58, 57, 53, 50, 48, 27, 47, 44, 1,
                             40, 2,  43, 0,  31, 45, 46, 29, 30, 42, 41, 3,  4,  39, 5,  38 } },
        /* 1NH */ { 12, 3, { 54, 23, 22, 19, 61, 17, 16, 63, 62, 60, 18, 20, 28, 49, 26, 51, 52, 25,
                             55, 58, 57, 53, 50, 48, 40, 2,  43, 0,  31, 45, 46, 29, 30, 42, 41, 3 } },
        /* 1NI */ { 14, 4, { 16, 63, 61, 17, 62, 19, 60, 18, 58, 20, 59, 22, 57, 21, 51, 50, 26, 49, 27,
                             48, 28, 47, 29, 46, 30, 45, 31, 44, 43, 0,  1,  42, 2,  41, 3,  40, 4,  39,
                             5,  38, 6,  37, 13, 34, 9,  35, 8,  36, 7,  10, 33, 11, 32, 12, 15, 14 } },
        /* 1NJ */ { 10, 4, { 16, 63, 61, 17, 62, 19, 60, 18, 58, 20, 51, 50, 26, 49, 27, 48, 28, 47, 29, 46,
                             43, 0,  1,  42, 2,  41, 3,  40, 4,  39, 13, 34, 9,  35, 8,  36, 7,  10, 33, 11 } },
        /* 1NK */ { 10, 4, { 16, 63, 61, 17, 62, 19, 60, 18, 58, 20, 51, 50, 26, 49, 27, 48, 28, 47, 29, 46,
                             43, 0,  1,  42, 2,  41, 3,  40, 4,  39, 13, 34, 9,  35, 8,  36, 7,  10, 33, 11 } },
        /* 1NL */ { 5, 18, { -1, -1, 63, 17, 62, -1, -1, 16, 19, 18, -1, -1, 61, 60, 59, -1, -1, 58, 20, 57, -1, -1, 22,
                             21, 23, -1, -1, 55, 56, 54, -1, -1, 25, 24, 53, -1, -1, 51, 50, 52, -1, -1, 49, 27, 26, -1,
                             48, 28, 29, 47, -1, 46, 31, 45, 30, -1, 43, 2,  42, 44, -1, 0,  3,  4,  1,  -1, 40, 6,  38,
                             41, -1, 37, 34, 36, 39, 5,  9,  33, 10, 7,  8,  11, 13, 14, 35, -1, -1, 12, 15, 32 } },
        /* 1NM */ { 5, 15, { -1, 61, -1, -1, -1, -1, 19, 16, 63, 62, -1, 22, 58, 17, 18, -1, 55, 57, 59,
                             60, -1, 54, 23, 21, 20, -1, 51, 25, 24, 56, -1, 49, 26, 52, 53, -1, 48, 28,
                             27, 50, 29, 46, 45, 30, 47, 31, 43, 2,  0,  44, 5,  39, 40, 42, 1,  8,  36,
                             37, 4,  41, 33, 10, 34, 6,  3,  12, 32, 11, 35, 38, 13, 14, 15, 9,  7 } },
        /* 1NN */ { 5, 14, { 61, 17, 16, 63, -1, 19, 60, 18, 62, -1, 21, 58, 20, 59, -1, 57, 23, 22,
                             56, -1, 55, 24, 54, 52, 53, 25, 51, 26, 50, 27, 49, 48, 28, 47, 29, 46,
                             45, 31, 30, 44, 0,  43, 42, 1,  41, 2,  3,  40, 4,  39, 5,  6,  37, 38,
                             7,  36, 8,  9,  35, 10, -1, 34, 33, 32, 12, -1, 11, 13, 15, 14 } },
        /* 1NG */ { 16, 1, { 13, 11, 6, 37, 8, 9, 34, 33, 12, 36, 7, 35, 10, 32, 14, 15 } },
        /* 1NH */ { 12, 1, { 13, 11, 6, 37, 8, 9, 34, 33, 12, 36, 7, 35 } } },
      /* PS */
      { { 0.63, 0.42 }, { 0.63, 0.84 }, { 0.63, 1.68 }, { 0.63, 3.36 } }
    };
  }
}
class SegmentationCreatorRegisterCreateSegType0
{
 public:
  SegmentationCreatorRegisterCreateSegType0() { registerSegmentationCreator(0, createSegType0); }
} aSegmentationCreatorRegisterCreateSegType0;

} // namespace impl2
} // namespace mapping
} // namespace mch
} // namespace o2