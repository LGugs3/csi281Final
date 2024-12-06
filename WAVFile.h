//
// Created by liamg on 12/4/2024.
//

#ifndef WAVHEADER_H
#define WAVHEADER_H

#include <bitset>
#include <cmath>
#include <complex>
#include <cstdint>
#include <valarray>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cassert>
#include <limits>

using Complex = std::complex<double>;
using CArray = std::valarray<Complex>;

class WAVFile
{
    private:
    struct WAVHeader
    {
        char chunkID[4];
        uint32_t chunkSize;
        char format[4];

        char subchunk1ID[4];
        uint32_t subchunk1Size;
        short audioFormat;
        short numChannels;
        uint32_t sampleRate;
        uint32_t byteRate;
        short blockAlign;
        short bitsPerSample;

        char subchunk2ID[4];
        uint32_t subchunk2Size;

    public:
        explicit WAVHeader(const std::string& file);
    };

    WAVHeader header;

    std::string fileName;
    std::vector<std::pair<int16_t, int16_t>> samples;
    std::vector<Complex> leftData;
    std::vector<Complex> rightData;

    public:
    explicit WAVFile(const std::string& fileName): header(fileName), fileName(fileName)
    {
        assert(header.numChannels <= 2);
        assert(header.numChannels > 0);
        assert(header.bitsPerSample == 16);
        samples = std::vector<std::pair<int16_t, int16_t>>();
        leftData = std::vector<Complex>();
        rightData = std::vector<Complex>();
    }
    void readData();
    void printData();
    void initiateFFT();
    void runBluesteinFFT(std::vector<Complex>& vec, bool inverse);
    std::vector<std::complex<double>> convolve(std::vector<Complex> xvec, std::vector<Complex> yvec);
    void transform(std::vector<Complex>& vec, bool inverse);
    void transformRadix2(std::vector<Complex> &vec, bool inverse);
    static size_t reverseBits(size_t val, int width);
    void createGraph();
};



#endif //WAVHEADER_H
