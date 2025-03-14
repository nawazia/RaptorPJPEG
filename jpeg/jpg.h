//
// Created by Ibrahim Nawaz on 07/01/2025.
//

#ifndef JPG_H
#define JPG_H
#include <cstring>
#include <map>

typedef unsigned char byte;
typedef unsigned int uint;

// Start/End of Image markers
const byte SOI = 0xD8;
const byte EOI = 0xD9;

// Start of Frame markers, non-differential, Huffman coding
const byte SOF0 = 0xC0; // Baseline DCT
const byte SOF1 = 0xC1; // Extended sequential DCT
const byte SOF2 = 0xC2; // Progressive DCT
const byte SOF3 = 0xC3; // Lossless (sequential)
// Start of Frame markers, differential, Huffman coding
const byte SOF5 = 0xC5; // Differential sequential DCT
const byte SOF6 = 0xC6; // Differential progressive DCT
const byte SOF7 = 0xC7; // Differential lossless (sequential)
// Start of Frame markers, non-differential, arithmetic coding
const byte SOF9 = 0xC9; // Extended sequential DCT
const byte SOF10 = 0xCA; // Progressive DCT
const byte SOF11 = 0xCB; // Lossless (sequential)
// Start of Frame markers, differential, arithmetic coding
const byte SOF13 = 0xCD; // Differential sequential DCT
const byte SOF14 = 0xCE; // Differential progressive DCT
const byte SOF15 = 0xCF; // Differential lossless (sequential)

// DRI, defined restart interval
const byte DRI = 0xDD;

// DHT, huffman tables
const byte DHT = 0xC4;

// Start of Scan
const byte SOS = 0xDA;

// APPN markers
const byte APP0 = 0xE0;
const byte APP1 = 0xE1;
const byte APP2 = 0xE2;
const byte APP3 = 0xE3;
const byte APP4 = 0xE4;
const byte APP5 = 0xE5;
const byte APP6 = 0xE6;
const byte APP7 = 0xE7;
const byte APP8 = 0xE8;
const byte APP9 = 0xE9;
const byte APP10 = 0xEA;
const byte APP11 = 0xEB;
const byte APP12 = 0xEC;
const byte APP13 = 0xED;
const byte APP14 = 0xEE;
const byte APP15 = 0xEF;

// DQT - can define multiple tables within the same marker, defined by marker length.
const byte DQT = 0xDB;

// Restart interval Markers
const byte RST0 = 0xD0;
const byte RST1 = 0xD1;
const byte RST2 = 0xD2;
const byte RST3 = 0xD3;
const byte RST4 = 0xD4;
const byte RST5 = 0xD5;
const byte RST6 = 0xD6;
const byte RST7 = 0xD7;

// Misc markers
const byte COM = 0xFE;
const byte TEM = 0x01;

const byte DNL = 0xDC; // Define Number of Lines
const byte DHP = 0xDE; // Define Hierarchical Progression
const byte EXP = 0xDF; // Expand Reference Component(s)
const byte DAC = 0xCC; // Define Arithmetic Coding Conditioning(s)

const byte JPG0 = 0xF0;
const byte JPG1 = 0xF1;
const byte JPG2 = 0xF2;
const byte JPG3 = 0xF3;
const byte JPG4 = 0xF4;
const byte JPG5 = 0xF5;
const byte JPG6 = 0xF6;
const byte JPG7 = 0xF7;
const byte JPG8 = 0xF8;
const byte JPG9 = 0xF9;
const byte JPG10 = 0xFA;
const byte JPG11 = 0xFB;
const byte JPG12 = 0xFC;
const byte JPG13 = 0xFD;


struct QuantizationTable {
    uint table[64] = { 0 };
    bool set = false;
};

struct HuffmanTable {
    byte offsets[17] = { 0 };
    byte symbols[162] = { 0 };
    uint codes[162] = { 0 };
    bool set = false;
};

struct ColorComponent {
    byte horizontalSamplingFactor = 1;
    byte verticalSamplingFactor = 1;
    byte quantizationTableID = 0;
    byte huffmanDCTableID = 0;
    byte huffmanACTableID = 0;
    // bool used = false;
    bool usedInFrame = false;
    bool usedInScan = false;
};

struct Block {
    union {
        int y[64] = { 0 };
        int r[64];
    };
    union {
        int cb[64] = { 0 };
        int g[64];
    };
    union {
        int cr[64] = { 0 };
        int b[64];
    };
    int* operator[](uint i) {
        switch (i) {
            case 0:
                return y;
            case 1:
                return cb;
            case 2:
                return cr;
            default:
                return nullptr;
        }
    }
};

struct JPGImage {
    QuantizationTable quantizationTables[4];
    HuffmanTable huffmanDCTables[4];
    HuffmanTable huffmanACTables[4];

    std::vector<std::int32_t> mapping;
    std::vector<std::vector<std::uint32_t>> validBits;

    byte frameType = 0;
    uint height = 0;
    uint width = 0;
    byte numComponents = 0;
    bool zeroBased = false;

    byte componentsInScan = 0;
    byte startOfSelection = 0;
    byte endOfSelection = 63;
    byte successiveApproximationHigh = 0;
    byte successiveApproximationLow = 0;

    uint restartInterval = 0;
    uint validityUpperBound = 0;
    int chunkSize;
    float avgBlockSize = 0;

    ColorComponent colorComponents[3];

    Block* blocks = nullptr;

    bool valid = true;
    bool validScan = true;

    uint blockHeight = 0;
    uint blockWidth = 0;
    uint blockHeightReal = 0;
    uint blockWidthReal = 0;

    // max sampling factor from all ColorComponent
    uint horizontalSamplingFactor = 1;
    uint verticalSamplingFactor = 1;
};

// IDCT scaling factors
const float m0 = 2.0 * std::cos(1.0 / 16.0 * 2.0 * M_PI);
const float m1 = 2.0 * std::cos(2.0 / 16.0 * 2.0 * M_PI);
const float m3 = 2.0 * std::cos(2.0 / 16.0 * 2.0 * M_PI);
const float m5 = 2.0 * std::cos(3.0 / 16.0 * 2.0 * M_PI);
const float m2 = m0 - m5;
const float m4 = m0 + m5;

const float s0 = std::cos(0.0 / 16.0 * M_PI) / std::sqrt(8);
const float s1 = std::cos(1.0 / 16.0 * M_PI) / 2.0;
const float s2 = std::cos(2.0 / 16.0 * M_PI) / 2.0;
const float s3 = std::cos(3.0 / 16.0 * M_PI) / 2.0;
const float s4 = std::cos(4.0 / 16.0 * M_PI) / 2.0;
const float s5 = std::cos(5.0 / 16.0 * M_PI) / 2.0;
const float s6 = std::cos(6.0 / 16.0 * M_PI) / 2.0;
const float s7 = std::cos(7.0 / 16.0 * M_PI) / 2.0;

const byte zigZagMap[] = {
    0,   1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63
};

const std::unordered_map<byte, int> RSTMap =
    {{0xD0, 0},
    {0xD1, 1},
    {0xD2, 2},
    {0xD3, 3},
    {0xD4, 4},
    {0xD5, 5},
    {0xD6, 6},
    {0xD7, 7}};

#endif //JPG_H
