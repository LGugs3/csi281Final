//
// Created by liamg on 11/22/2024.
//
#include <cmath>
#include <complex>
#include <valarray>
#include <iostream>
#include "AudioFile.h"

using Complex = std::complex<double>;
using CArray = std::valarray<Complex>;

void fft(CArray& arr)
{
    const size_t N = arr.size();
    if (N <= 1) return;

    //divide
    CArray even = arr[std::slice(0,N / 2, 2)];
    CArray odd = arr[std::slice(1,N / 2, 2)];

    fft(even);
    fft(odd);

    //combine
    for(size_t i = 0; i < N / 2; i++)
    {
        Complex t = std::polar(1.0, -2 * M_PI * i / N) * odd[i];
        arr[i] = even[i] + t;
        arr[i + N / 2] = even[i] - t;
    }
}

int main()
{
    AudioFile<double> audioFile;
    audioFile.load("basicScale.wav");

    audioFile.printSummary();

    // Complex test[8];
    // for(int i = 0; i < 4; i++)
    // {
    //     test[i] = 1.0;
    //     test[i + 4] = 0.0;
    // }
    // CArray data(test, 8);
    //
    // for(int i = 0; i < 8; i++)
    // {
    //     std::cout << data[i] << std::endl;
    // }
    // std::cout << std::endl << std::endl;
    // fft(data);
    //
    // for (int i = 0 ; i < 8; i++)
    // {
    //     std::cout << data[i] << std::endl;
    // }
    return 0;
}