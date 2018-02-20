//
// This file has been generated. Do not modify it by hand or your changes might be lost.
//
#include "SegmentationCreator.h"

namespace o2
{
namespace mch
{
namespace mapping
{
namespace impl3
{
Segmentation* createSegType3(bool isBendingPlane)
{
  if (isBendingPlane) {
    return new Segmentation{
      3,
      true,
      /* PG */
      { { 4, 20, 0, 40, -20 },     { 5, 21, 0, 45, -20 },     { 6, 21, 0, 50, -20 },    { 7, 11, 0, 55, -17.5 },
        { 101, 13, 1, -75, -20 },  { 102, 14, 1, -70, -20 },  { 103, 19, 1, -65, -20 }, { 104, 26, 1, -60, -20 },
        { 105, 18, 1, -50, -20 },  { 109, 17, 1, -40, -20 },  { 110, 25, 1, -35, -20 }, { 111, 19, 1, -25, -20 },
        { 112, 26, 1, -20, -20 },  { 113, 18, 1, -10, -20 },  { 118, 17, 0, 0, -20 },   { 119, 25, 0, 2.5, -20 },
        { 120, 19, 0, 7.5, -20 },  { 121, 26, 0, 10, -20 },   { 122, 18, 0, 15, -20 },  { 123, 17, 0, 20, -20 },
        { 124, 25, 0, 22.5, -20 }, { 125, 19, 0, 27.5, -20 }, { 126, 26, 0, 30, -20 },  { 127, 18, 0, 35, -20 },
        { 204, 16, 1, -50, 0 },    { 205, 24, 1, -60, 0 },    { 206, 22, 1, -65, 4 },   { 207, 27, 1, -75, 0 },
        { 211, 16, 1, -10, 0 },    { 212, 24, 1, -20, 0 },    { 213, 22, 1, -25, 4 },   { 214, 23, 1, -35, 0 },
        { 215, 15, 1, -40, 0 },    { 223, 16, 0, 35, 0 },     { 224, 24, 0, 30, 0 },    { 225, 22, 0, 27.5, 4 },
        { 226, 23, 0, 22.5, 0 },   { 227, 15, 0, 20, 0 },     { 228, 16, 0, 15, 0 },    { 229, 24, 0, 10, 0 },
        { 230, 22, 0, 7.5, 4 },    { 231, 23, 0, 2.5, 0 },    { 232, 15, 0, 0, 0 },     { 401, 12, 0, 75, -7 },
        { 402, 0, 0, 72.5, -7.5 }, { 403, 1, 0, 70, -8 },     { 404, 2, 0, 67.5, 1 },   { 405, 3, 0, 65, -8.5 },
        { 406, 4, 0, 62.5, -10 },  { 407, 5, 0, 60, -11 },    { 408, 6, 0, 55, -4 },    { 409, 7, 0, 52.5, -4 },
        { 410, 8, 0, 50, -4 },     { 411, 10, 0, 45, -4 },    { 412, 7, 0, 42.5, -4 },  { 413, 9, 0, 40, -4 } },
      /* PGT */
      { /* A10 */ { 2, 55, { 20, -1, 19, -1, 18, -1, 17, -1, 16, -1, 15, -1, 14, -1, 13, -1, 12, -1, 11, -1, 10, -1,
                             9,  -1, 8,  47, 7,  46, 6,  45, 5,  44, 4,  43, 3,  42, 2,  41, 1,  40, 0,  39, -1, 38,
                             -1, 37, -1, 36, -1, 35, -1, 34, -1, 33, -1, 32, -1, 63, -1, 62, -1, 61, -1, 60, -1, 59,
                             -1, 58, -1, 57, -1, 56, -1, 55, -1, 54, -1, 53, -1, 52, -1, 51, -1, 50, -1, 49, -1, 48,
                             -1, 21, -1, 22, -1, 23, -1, 24, -1, 25, -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31 } },
        /* A11 */ { 2, 56, { 29, -1, 28, -1, 27, -1, 26, -1, 25, -1, 24, -1, 23, -1, 22, -1, 21, -1, 20, -1, 19, -1, 18,
                             -1, 17, -1, 16, -1, 15, -1, 14, -1, 13, -1, 12, -1, 11, -1, 10, -1, 9,  -1, 8,  -1, 7,  47,
                             6,  46, 5,  45, 4,  44, 3,  43, 2,  42, 1,  41, 0,  40, -1, 39, -1, 38, -1, 37, -1, 36, -1,
                             35, -1, 34, -1, 33, -1, 32, -1, 63, -1, 62, -1, 61, -1, 60, -1, 59, -1, 58, -1, 57, -1, 56,
                             -1, 55, -1, 54, -1, 53, -1, 52, -1, 51, -1, 50, -1, 49, -1, 48, -1, 30, -1, 31 } },
        /* A12 */ { 2, 38, { 47, -1, 46, -1, 45, -1, 44, -1, 43, -1, 42, -1, 31, -1, 30, -1, 29, -1, 28,
                             -1, 27, -1, 26, -1, 25, 41, 24, 40, 23, 39, 22, 38, 21, 37, 20, 36, 19, 35,
                             18, 34, 17, 33, 16, 32, 15, 63, 14, 62, 13, 61, 12, 60, 11, 59, 10, 58, 9,
                             57, 8,  56, 7,  55, 6,  54, 5,  53, 4,  52, 3,  51, 2,  50, 1,  49, 0,  48 } },
        /* A13 */ { 2, 57, { -1, 61, -1, 62, -1, 63, -1, 32, -1, 33, -1, 34, -1, 35, -1, 36, -1, 37, -1, 38, -1, 39, -1,
                             40, 60, 41, 59, 42, 58, 43, 57, 44, 56, 45, 55, 46, 54, 47, 53, -1, 52, -1, 51, -1, 50, -1,
                             49, -1, 48, -1, 31, -1, 30, -1, 29, -1, 28, -1, 27, -1, 26, -1, 25, -1, 24, -1, 23, -1, 22,
                             -1, 21, -1, 20, -1, 19, -1, 18, -1, 17, -1, 16, -1, 15, -1, 14, -1, 13, -1, 12, -1, 11, -1,
                             10, -1, 9,  -1, 8,  -1, 7,  -1, 6,  -1, 5,  -1, 4,  -1, 3,  -1, 2,  -1, 1,  -1, 0,  -1 } },
        /* A14 */ { 2, 60, { -1, 17, -1, 18, -1, 19, -1, 20, -1, 21, -1, 22, -1, 23, -1, 24, -1, 25, -1, 26,
                             -1, 27, 16, 28, 15, 29, 14, 30, 13, 31, 12, -1, 11, -1, 10, -1, 9,  -1, 8,  -1,
                             7,  -1, 6,  -1, 5,  -1, 4,  -1, 3,  -1, 2,  -1, 1,  -1, 0,  -1, 47, -1, 46, -1,
                             45, -1, 44, -1, 43, -1, 42, -1, 41, -1, 40, -1, 39, -1, 38, -1, 37, -1, 36, -1,
                             35, -1, 34, -1, 33, -1, 32, -1, 63, -1, 62, -1, 61, -1, 60, -1, 59, -1, 58, -1,
                             57, -1, 56, -1, 55, -1, 54, -1, 53, -1, 52, -1, 51, -1, 50, -1, 49, -1, 48, -1 } },
        /* A15 */ { 2, 62, { -1, 35, -1, 36, -1, 37, -1, 38, -1, 39, -1, 40, -1, 41, -1, 42, -1, 43, -1, 44, -1,
                             45, 34, 46, 33, 47, 32, -1, 63, -1, 62, -1, 61, -1, 60, -1, 59, -1, 58, -1, 57, -1,
                             56, -1, 55, -1, 54, -1, 53, -1, 52, -1, 51, -1, 50, -1, 49, -1, 48, -1, 0,  -1, 1,
                             -1, 2,  -1, 3,  -1, 4,  -1, 5,  -1, 6,  -1, 7,  -1, 8,  -1, 9,  -1, 10, -1, 11, -1,
                             12, -1, 13, -1, 14, -1, 15, -1, 16, -1, 17, -1, 18, -1, 19, -1, 20, -1, 21, -1, 22,
                             -1, 23, -1, 24, -1, 25, -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31, -1 } },
        /* A16 */ { 2, 48, { -1, 16, -1, 17, -1, 18, -1, 19, -1, 20, -1, 21, -1, 22, -1, 23, -1, 24, -1, 25,
                             -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31, -1, 47, -1, 46, -1, 45, -1, 44,
                             -1, 43, -1, 42, -1, 41, -1, 40, -1, 39, -1, 38, -1, 37, -1, 36, -1, 35, -1, 34,
                             -1, 33, -1, 32, 0,  63, 1,  62, 2,  61, 3,  60, 4,  59, 5,  58, 6,  57, 7,  56,
                             8,  55, 9,  54, 10, 53, 11, 52, 12, 51, 13, 50, 14, 49, 15, 48 } },
        /* A17 */ { 2, 48, { -1, 48, -1, 49, -1, 50, -1, 51, -1, 52, -1, 53, -1, 54, -1, 55, -1, 56, -1, 57,
                             -1, 58, -1, 59, -1, 60, -1, 61, -1, 62, -1, 63, 0,  32, 1,  33, 2,  34, 3,  35,
                             4,  36, 5,  37, 6,  38, 7,  39, 8,  40, 9,  41, 10, 42, 11, 43, 12, 44, 13, 45,
                             14, 46, 15, 47, 16, -1, 17, -1, 18, -1, 19, -1, 20, -1, 21, -1, 22, -1, 23, -1,
                             24, -1, 25, -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31, -1 } },
        /* A18 */ { 2, 48, { 63, 32, 62, 33, 61, 34, 60, 35, 59, 36, 58, 37, 57, 38, 56, 39, 55, 40, 54, 41,
                             53, 42, 52, 43, 51, 44, 50, 45, 49, 46, 48, 47, 0,  -1, 1,  -1, 2,  -1, 3,  -1,
                             4,  -1, 5,  -1, 6,  -1, 7,  -1, 8,  -1, 9,  -1, 10, -1, 11, -1, 12, -1, 13, -1,
                             14, -1, 15, -1, 16, -1, 17, -1, 18, -1, 19, -1, 20, -1, 21, -1, 22, -1, 23, -1,
                             24, -1, 25, -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31, -1 } },
        /* A19 */ { 2, 48, { 63, 32, 62, 33, 61, 34, 60, 35, 59, 36, 58, 37, 57, 38, 56, 39, 55, 40, 54, 41,
                             53, 42, 52, 43, 51, 44, 50, 45, 49, 46, 48, 47, 31, -1, 30, -1, 29, -1, 28, -1,
                             27, -1, 26, -1, 25, -1, 24, -1, 23, -1, 22, -1, 21, -1, 20, -1, 19, -1, 18, -1,
                             17, -1, 16, -1, 15, -1, 14, -1, 13, -1, 12, -1, 11, -1, 10, -1, 9,  -1, 8,  -1,
                             7,  -1, 6,  -1, 5,  -1, 4,  -1, 3,  -1, 2,  -1, 1,  -1, 0,  -1 } },
        /* A20 */ { 2, 48, { -1, 16, -1, 17, -1, 18, -1, 19, -1, 20, -1, 21, -1, 22, -1, 23, -1, 24, -1, 25,
                             -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31, -1, 48, -1, 49, -1, 50, -1, 51,
                             -1, 52, -1, 53, -1, 54, -1, 55, -1, 56, -1, 57, -1, 58, -1, 59, -1, 60, -1, 61,
                             -1, 62, -1, 63, 0,  32, 1,  33, 2,  34, 3,  35, 4,  36, 5,  37, 6,  38, 7,  39,
                             8,  40, 9,  41, 10, 42, 11, 43, 12, 44, 13, 45, 14, 46, 15, 47 } },
        /* A8 */ { 3, 27, { 48, -1, -1, 49, -1, -1, 50, -1, -1, 51, -1, -1, 52, -1, -1, 53, 0,  -1, 54, 1,  -1,
                            55, 2,  -1, 56, 3,  -1, 57, 5,  4,  58, 7,  6,  59, 9,  8,  60, 11, 10, 61, 13, 12,
                            62, 15, 14, 63, 17, 16, 32, 19, 18, 33, 21, 20, 34, 23, 22, 35, 25, 24, 36, 27, 26,
                            37, 29, 28, 39, 38, 30, 41, 40, 31, 43, 42, -1, 45, 44, -1, 47, 46, -1 } },
        /* A9 */ { 2, 54, { 42, -1, 40, 47, 38, 46, 36, 45, 34, 44, 32, 43, 62, 41, 60, 39, 58, 37, 56, 35, 54, 33,
                            -1, 63, -1, 61, -1, 59, -1, 57, -1, 55, -1, 53, -1, 52, -1, 51, -1, 50, -1, 49, -1, 48,
                            -1, 31, -1, 30, -1, 29, -1, 28, -1, 27, -1, 26, -1, 25, -1, 24, -1, 23, -1, 22, -1, 21,
                            -1, 20, -1, 19, -1, 18, -1, 17, -1, 16, -1, 15, -1, 14, -1, 13, -1, 12, -1, 11, -1, 10,
                            -1, 9,  -1, 8,  -1, 7,  -1, 6,  -1, 5,  -1, 4,  -1, 3,  -1, 2,  -1, 1,  -1, 0 } },
        /* I1 */ { 1, 64, { 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 63, 62, 61, 60, 59, 58,
                            57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20,
                            19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9,  8,  7,  6,  5,  4,  3,  2,  1,  0 } },
        /* L18 */ { 2, 40, { 63, -1, 62, -1, 61, -1, 60, -1, 59, -1, 58, -1, 57, -1, 56, -1, 55, -1, 54, -1,
                             53, -1, 52, -1, 51, -1, 50, -1, 49, -1, 48, -1, 31, -1, 30, -1, 29, -1, 28, -1,
                             27, -1, 26, -1, 25, -1, 24, -1, 23, -1, 22, -1, 21, -1, 20, -1, 19, -1, 18, -1,
                             17, -1, 16, -1, 15, 0,  14, 1,  13, 2,  12, 3,  11, 4,  10, 5,  9,  6,  8,  7 } },
        /* L5 */ { 2, 40, { 55, 56, 54, 57, 53, 58, 52, 59, 51, 60, 50, 61, 49, 62, 48, 63, 31, 32, 30, 33,
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
        /* O23 */ { 2, 32, { 47, 0,  46, 1,  45, 2,  44, 3,  43, 4,  42, 5,  41, 6,  40, 7,  39, 8,  38, 9,  37, 10,
                             36, 11, 35, 12, 34, 13, 33, 14, 32, 15, 63, 16, 62, 17, 61, 18, 60, 19, 59, 20, 58, 21,
                             57, 22, 56, 23, 55, 24, 54, 25, 53, 26, 52, 27, 51, 28, 50, 29, 49, 30, 48, 31 } },
        /* O24 */ { 2, 32, { 48, 0,  49, 1,  50, 2,  51, 3,  52, 4,  53, 5,  54, 6,  55, 7,  56, 8,  57, 9,  58, 10,
                             59, 11, 60, 12, 61, 13, 62, 14, 63, 15, 32, 16, 33, 17, 34, 18, 35, 19, 36, 20, 37, 21,
                             38, 22, 39, 23, 40, 24, 41, 25, 42, 26, 43, 27, 44, 28, 45, 29, 46, 30, 47, 31 } },
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
                            33, -1, 45, 34, -1, 44, 35, -1, 43, 36, -1, 42, 37, -1, 41, 38, -1, 40, 39, -1 } },
        /* Z5 */ { 3, 40, { -1, 39, 40, -1, 38, 41, -1, 37, 42, -1, 36, 43, -1, 35, 44, -1, 34, 45, -1, 33,
                            46, -1, 32, 47, -1, 63, -1, -1, 62, -1, -1, 61, -1, -1, 60, -1, -1, 59, -1, -1,
                            58, -1, -1, 57, -1, -1, 56, -1, -1, 55, -1, -1, 54, -1, -1, 53, -1, -1, 52, -1,
                            -1, 51, -1, -1, 50, -1, -1, 49, -1, -1, 48, -1, 0,  31, -1, 1,  30, -1, 2,  29,
                            -1, 3,  28, -1, 4,  27, -1, 5,  26, -1, 6,  25, -1, 7,  24, -1, 8,  23, -1, 9,
                            22, -1, 10, 21, -1, 11, 20, -1, 12, 19, -1, 13, 18, -1, 14, 17, -1, 15, 16, -1 } } },
      /* PS */
      { { 2.5, 0.5 }, { 5, 0.5 } }
    };
  } else {
    return new Segmentation{
      3,
      false,
      /* PG */
      { { 1025, 0, 0, 51.42856979, -20 },
        { 1026, 11, 0, 45.7142868, -20 },
        { 1027, 11, 0, 40, -20 },
        { 1130, 13, 1, -51.42856979, -20 },
        { 1131, 13, 1, -62.8571434, -20 },
        { 1132, 13, 1, -74.2857132, -20 },
        { 1138, 16, 1, -10, -20 },
        { 1139, 18, 1, -20, -20 },
        { 1140, 17, 1, -31.4285717, -20 },
        { 1141, 15, 1, -40, -20 },
        { 1152, 10, 0, 34.2857132, -20 },
        { 1153, 10, 0, 28.5714283, -20 },
        { 1154, 10, 0, 22.8571434, -20 },
        { 1155, 10, 0, 17.1428566, -20 },
        { 1156, 10, 0, 11.4285717, -20 },
        { 1157, 10, 0, 5.714285851, -20 },
        { 1158, 10, 0, -3.996802889e-15, -20 },
        { 1225, 14, 1, -74.2857132, 0 },
        { 1226, 14, 1, -62.8571434, 0 },
        { 1227, 14, 1, -51.42856979, 0 },
        { 1232, 7, 1, -40, 0 },
        { 1233, 14, 1, -25.7142849, 0 },
        { 1234, 8, 1, -14.28571415, 0 },
        { 1240, 9, 0, -3.996802889e-15, 0 },
        { 1241, 9, 0, 5.714285851, 0 },
        { 1242, 9, 0, 11.4285717, 0 },
        { 1243, 9, 0, 17.1428566, 0 },
        { 1244, 9, 0, 22.8571434, 0 },
        { 1245, 9, 0, 28.5714283, 0 },
        { 1246, 9, 0, 34.2857132, 0 },
        { 1325, 12, 0, 40, 0 },
        { 1326, 12, 0, 45.7142868, 0 },
        { 1327, 12, 0, 51.42856979, 0 },
        { 1328, 1, 0, 57.1428566, -15 },
        { 1329, 2, 0, 60, -12.5 },
        { 1330, 3, 0, 63.57143021, -10 },
        { 1331, 4, 0, 67.14286041, -10 },
        { 1332, 5, 0, 71.42857361, -7.5 },
        { 1333, 6, 0, 75.7142868, -7.5 } },
      /* PGT */
      { /* A1 */ {
          9, 8, { 7,  15, 23, 31, 55, 63, -1, -1, -1, 6,  14, 22, 30, 54, 62, 38, 45, 47, 5,  13, 21, 29, 53, 61,
                  37, 44, 46, 4,  12, 20, 28, 52, 60, 36, 43, -1, 3,  11, 19, 27, 51, 59, 35, 42, -1, 2,  10, 18,
                  26, 50, 58, 34, 41, -1, 1,  9,  17, 25, 49, 57, 33, 40, -1, 0,  8,  16, 24, 48, 56, 32, 39, -1 } },
        /* A2 */ { 5, 14, { -1, 34, 52, 22, 8,  47, 33, 51, 21, 7,  46, 32, 50, 20, 6,  45, 63, 49,
                            19, 5,  44, 62, 48, 18, 4,  43, 61, 31, 17, 3,  42, 60, 30, 16, 2,  41,
                            59, 29, 15, 1,  40, 58, 28, 14, 0,  39, 57, 27, 13, -1, 38, 56, 26, 12,
                            -1, 37, 55, 25, 11, -1, 36, 54, 24, 10, -1, 35, 53, 23, 9,  -1 } },
        /* A3 */ { 6, 13, { -1, 42, 61, 48, 19, 6,  -1, 41, 60, 31, 18, 5,  -1, 40, 59, 30, 17, 4,  -1, 39,
                            58, 29, 16, 3,  -1, 38, 57, 28, 15, 2,  -1, 37, 56, 27, 14, 1,  -1, 36, 55, 26,
                            13, 0,  -1, 35, 54, 25, 12, -1, 47, 34, 53, 24, 11, -1, 46, 33, 52, 23, 10, -1,
                            45, 32, 51, 22, 9,  -1, 44, 63, 50, 21, 8,  -1, 43, 62, 49, 20, 7,  -1 } },
        /* A4 */
        { 6, 12, { -1, 41, 61, 49, 21, 9, -1, 40, 60, 48, 20, 8, -1, 39, 59, 31, 19, 7,  -1, 38, 58, 30, 18, 6,
                   -1, 37, 57, 29, 17, 5, -1, 36, 56, 28, 16, 4, 47, 35, 55, 27, 15, 3,  46, 34, 54, 26, 14, 2,
                   45, 33, 53, 25, 13, 1, 44, 32, 52, 24, 12, 0, 43, 63, 51, 23, 11, -1, 42, 62, 50, 22, 10, -1 } },
        /* A5 */ { 7, 12, { -1, 45, 33, 53, -1, -1, -1, -1, 44, 32, 52, 25, 14, 3,  -1, 43, 63, 51, 24, 13, 2,
                            -1, 42, 62, 50, 23, 12, 1,  -1, 41, 61, 49, 22, 11, 0,  -1, 40, 60, 48, 21, 10, -1,
                            -1, 39, 59, 31, 20, 9,  -1, -1, 38, 58, 30, 19, 8,  -1, -1, 37, 57, 29, 18, 7,  -1,
                            -1, 36, 56, 28, 17, 6,  -1, 47, 35, 55, 27, 16, 5,  -1, 46, 34, 54, 26, 15, 4,  -1 } },
        /* A6 */ { 7, 11, { -1, 40, 61, 50, 23, 12, 1,  -1, 39, 60, 49, 22, 11, 0,  -1, 38, 59, 48, 21, 10,
                            -1, -1, 37, 58, 31, 20, 9,  -1, 47, 36, 57, 30, 19, 8,  -1, 46, 35, 56, 29, 18,
                            7,  -1, 45, 34, 55, 28, 17, 6,  -1, 44, 33, 54, 27, 16, 5,  -1, 43, 32, 53, 26,
                            15, 4,  -1, 42, 63, 52, 25, 14, 3,  -1, 41, 62, 51, 24, 13, 2,  -1 } },
        /* A7 */ { 6, 11, { -1, 38, 59, 48, 21, 10, -1, 37, 58, 31, 20, 9,  47, 36, 57, 30, 19, 8,  46, 35, 56, 29,
                            18, 7,  45, 34, 55, 28, 17, 6,  44, 33, 54, 27, 16, 5,  43, 32, 53, 26, 15, 4,  42, 63,
                            52, 25, 14, 3,  41, 62, 51, 24, 13, 2,  40, 61, 50, 23, 12, 1,  39, 60, 49, 22, 11, 0 } },
        /* L3 */ { 20, 4, { 44, 40, 36, 32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                            45, 41, 37, 33, 61, 58, 55, 52, 49, 30, 27, 24, 21, 18, 15, 12, 9,  6,  3,  0,
                            46, 42, 38, 34, 62, 59, 56, 53, 50, 31, 28, 25, 22, 19, 16, 13, 10, 7,  4,  1,
                            47, 43, 39, 35, 63, 60, 57, 54, 51, 48, 29, 26, 23, 20, 17, 14, 11, 8,  5,  2 } },
        /* L4 */ { 20, 4, { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 12, 8,  4, 0,
                            45, 42, 39, 36, 33, 62, 59, 56, 53, 50, 31, 28, 25, 22, 19, 16, 13, 9,  5, 1,
                            46, 43, 40, 37, 34, 63, 60, 57, 54, 51, 48, 29, 26, 23, 20, 17, 14, 10, 6, 2,
                            47, 44, 41, 38, 35, 32, 61, 58, 55, 52, 49, 30, 27, 24, 21, 18, 15, 11, 7, 3 } },
        /* O1 */ { 8, 8, { 40, 32, 56, 48, 24, 16, 8,  0,  41, 33, 57, 49, 25, 17, 9,  1,  42, 34, 58, 50, 26, 18,
                           10, 2,  43, 35, 59, 51, 27, 19, 11, 3,  44, 36, 60, 52, 28, 20, 12, 4,  45, 37, 61, 53,
                           29, 21, 13, 5,  46, 38, 62, 54, 30, 22, 14, 6,  47, 39, 63, 55, 31, 23, 15, 7 } },
        /* O2 */ { 8, 8, { 7,  15, 23, 31, 55, 63, 39, 47, 6,  14, 22, 30, 54, 62, 38, 46, 5,  13, 21, 29, 53, 61,
                           37, 45, 4,  12, 20, 28, 52, 60, 36, 44, 3,  11, 19, 27, 51, 59, 35, 43, 2,  10, 18, 26,
                           50, 58, 34, 42, 1,  9,  17, 25, 49, 57, 33, 41, 0,  8,  16, 24, 48, 56, 32, 40 } },
        /* O21 */ { 8, 8, { 7,  15, 23, 31, 55, 63, 39, 47, 6,  14, 22, 30, 54, 62, 38, 46, 5,  13, 21, 29, 53, 61,
                            37, 45, 4,  12, 20, 28, 52, 60, 36, 44, 3,  11, 19, 27, 51, 59, 35, 43, 2,  10, 18, 26,
                            50, 58, 34, 42, 1,  9,  17, 25, 49, 57, 33, 41, 0,  8,  16, 24, 48, 56, 32, 40 } },
        /* O22 */ { 8, 8, { 47, 39, 63, 55, 31, 23, 15, 7,  46, 38, 62, 54, 30, 22, 14, 6,  45, 37, 61, 53, 29, 21,
                            13, 5,  44, 36, 60, 52, 28, 20, 12, 4,  43, 35, 59, 51, 27, 19, 11, 3,  42, 34, 58, 50,
                            26, 18, 10, 2,  41, 33, 57, 49, 25, 17, 9,  1,  40, 32, 56, 48, 24, 16, 8,  0 } },
        /* O3 */ { 16, 4, { 3,  7,  11, 15, 19, 23, 27, 31, 51, 55, 59, 63, 35, 39, 43, 47, 2,  6,  10, 14, 18, 22,
                            26, 30, 50, 54, 58, 62, 34, 38, 42, 46, 1,  5,  9,  13, 17, 21, 25, 29, 49, 53, 57, 61,
                            33, 37, 41, 45, 0,  4,  8,  12, 16, 20, 24, 28, 48, 52, 56, 60, 32, 36, 40, 44 } },
        /* O4 */ { 16, 4, { 44, 40, 36, 32, 60, 56, 52, 48, 28, 24, 20, 16, 12, 8,  4,  0,  45, 41, 37, 33, 61, 57,
                            53, 49, 29, 25, 21, 17, 13, 9,  5,  1,  46, 42, 38, 34, 62, 58, 54, 50, 30, 26, 22, 18,
                            14, 10, 6,  2,  47, 43, 39, 35, 63, 59, 55, 51, 31, 27, 23, 19, 15, 11, 7,  3 } },
        /* P3 */ { 14, 5, { 3,  7,  11, 15, 20, 25, 30, 51, 56, 61, 34, 39, 43, 47, 2,  6,  10, 14,
                            19, 24, 29, 50, 55, 60, 33, 38, 42, 46, 1,  5,  9,  13, 18, 23, 28, 49,
                            54, 59, 32, 37, 41, 45, 0,  4,  8,  12, 17, 22, 27, 48, 53, 58, 63, 36,
                            40, 44, -1, -1, -1, -1, 16, 21, 26, 31, 52, 57, 62, 35, -1, -1 } },
        /* P4 */ { 14, 5, { 3,  7,  12, 17, 22, 27, 48, 53, 58, 63, 35, 39, 43, 47, 2,  6,  11, 16,
                            21, 26, 31, 52, 57, 62, 34, 38, 42, 46, 1,  5,  10, 15, 20, 25, 30, 51,
                            56, 61, 33, 37, 41, 45, 0,  4,  9,  14, 19, 24, 29, 50, 55, 60, 32, 36,
                            40, 44, -1, -1, 8,  13, 18, 23, 28, 49, 54, 59, -1, -1, -1, -1 } },
        /* Q3 */ { 16, 5, { -1, -1, 6,  11, 16, 21, 26, 31, 51, 55, 59, 63, 35, 39, 43, 47, -1, -1, 5,  10,
                            15, 20, 25, 30, 50, 54, 58, 62, 34, 38, 42, 46, -1, -1, 4,  9,  14, 19, 24, 29,
                            49, 53, 57, 61, 33, 37, 41, 45, -1, -1, 3,  8,  13, 18, 23, 28, 48, 52, 56, 60,
                            32, 36, 40, 44, 0,  1,  2,  7,  12, 17, 22, 27, -1, -1, -1, -1, -1, -1, -1, -1 } },
        /* Q4 */ { 16, 5, { 3,  7,  11, 15, 19, 23, 27, 31, 52, 57, 62, 35, 40, 45, -1, -1, 2,  6,  10, 14,
                            18, 22, 26, 30, 51, 56, 61, 34, 39, 44, -1, -1, 1,  5,  9,  13, 17, 21, 25, 29,
                            50, 55, 60, 33, 38, 43, -1, -1, 0,  4,  8,  12, 16, 20, 24, 28, 49, 54, 59, 32,
                            37, 42, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 48, 53, 58, 63, 36, 41, 46, 47 } } },
      /* PS */
      { { 0.714285714, 2.5 }, { 0.714285714, 5 } }
    };
  }
}
class SegmentationCreatorRegisterCreateSegType3
{
 public:
  SegmentationCreatorRegisterCreateSegType3() { registerSegmentationCreator(3, createSegType3); }
} aSegmentationCreatorRegisterCreateSegType3;

} // namespace impl3
} // namespace mapping
} // namespace mch
} // namespace o2
