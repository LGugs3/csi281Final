//
// Created by liamg on 11/22/2024.
//

#include "WAVFile.h"

int main()
{

    WAVFile file("basicScale.wav");

    file.readData();
    std::cout << "Data read" << std::endl;
    file.initiateFFT();
    std::cout << "FFT finished" << std::endl;
    //file.printData();
    file.createGraph();
    std::cout << "Graph created" << std::endl;
    return 0;
}