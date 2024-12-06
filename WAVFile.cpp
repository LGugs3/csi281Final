//
// Created by liamg on 12/4/2024.
//
#include "WAVFile.h"

WAVFile::WAVHeader::WAVHeader(const std::string& file)
{
    std::ifstream fin(file, std::ios::binary);

    if (!fin.is_open())
    {
        std::cerr << "Can't open file" << std::endl;
    }

    fin.read(reinterpret_cast<char*>(this), sizeof(WAVFile));//forces cast binary to char* and read the first sizeof(WAVHeader) chars into header
    fin.close();

    if (strncmp(chunkID, "RIFF", 4) != 0 || strncmp(format, "WAVE", 4) != 0)//need to use stncmp b/c of char*
    {
        std::cerr << "Invalid format" << std::endl;
        exit(1);
    }
}

void WAVFile::readData()
{
    std::ifstream fin(fileName, std::ios::binary);
    if (!fin.is_open())
    {
        std::cerr << "Can't open file" << std::endl;
        return;
    }

    fin.seekg(44, std::ios::beg);

    std::pair<int16_t, int16_t> sample;
    while(!fin.eof())
    {
        fin.read(reinterpret_cast<char*>(&sample.first), sizeof(sample.first));
        fin.read(reinterpret_cast<char*>(&sample.second), sizeof(sample.second));

        samples.push_back(sample);
    }
}

void WAVFile::printData()
{
    std::cout <<"Left" << std::endl;
    for(auto it : leftData)
        std::cout << it.real() << std::endl;

    if(header.numChannels == 2)
    {
        std::cout <<"Right" << std::endl;
        for(auto it : rightData)
            std::cout << it.real() << std::endl;
    }
}

void WAVFile::initiateFFT()
{
    for(auto it : samples)
    {
        Complex var(it.first, 0);
        leftData.push_back(var);
    }
    runBluesteinFFT(leftData, false);

    if (header.numChannels > 1)
    {
        for(auto it : samples)
        {
            Complex var(it.second, 0);
            rightData.push_back(var);
        }
        runBluesteinFFT(rightData, false);
    }
}


void WAVFile::runBluesteinFFT(std::vector<Complex>& vec, bool inverse)
{
    size_t size = vec.size();
    size_t m = 1;
    //m needs to be a power of two such that m >= n * 2 + 1
    while (m / 2 <= size)
    {
        if (m > vec.max_size() / 2)
            throw std::length_error("Vector too long");

        m *= 2;
    }

    //trig tables
    std::vector<Complex> expTable(size);
    for(size_t i = 0; i < size; i++)
    {
        std::uintmax_t temp = i * i;
        temp %= size * 2;
        double angle = (inverse ? M_PI : -M_PI) * temp / size;
        expTable[i] = std::polar(1.0, angle);
    }

    //preprocessing
    std::vector<Complex> avec(m);
    for(size_t i = 0; i < size; i++)
        avec[i] = vec[i] * expTable[i];

    std::vector<Complex> bvec(m);
    bvec[0] = expTable[0];
    for(size_t i = 1; i < size; i++)
        bvec[i] = bvec[m - 1] = std::conj(expTable[i]);

    //Convolution
    std::vector<Complex> cvec = convolve(std::move(avec), std::move(bvec));

    //Post-processing
    for(size_t i = 0; i < size; i++)
        vec[i] = cvec[i] * expTable[i];
}

std::vector<std::complex<double>> WAVFile::convolve(std::vector<Complex> xvec, std::vector<Complex> yvec)
{
    size_t size = xvec.size();
    if (size != yvec.size())
    {
        std::cerr << "FFT convolve size mismatch" << std::endl;
    }

    transform(xvec, false);
    transform(yvec, false);
    for(size_t i = 0; i < size; i++)
        xvec[i] *= yvec[i];
    transform(xvec, true);
    for(size_t i = 0; i < size; i++)
        xvec[i] /= static_cast<double>(size);//scaling

    return xvec;

}

void WAVFile::transform(std::vector<Complex>& vec, bool inverse)
{
    size_t size = vec.size();
    if (size == 0) return;

    if ((size & (size - 1)) == 0)
        transformRadix2(vec, inverse);
    else
        runBluesteinFFT(vec, inverse);
}

void WAVFile::transformRadix2(std::vector<Complex> &vec, bool inverse) {
    // Length variables
    size_t n = vec.size();
    int levels = 0;  // Compute levels = floor(log2(n))
    for (size_t temp = n; temp > 1U; temp >>= 1)
        levels++;
    if (static_cast<size_t>(1U) << levels != n)
        throw std::domain_error("Length is not a power of 2");

    // Trigonometric table
    std::vector<Complex> expTable(n / 2);
    for (size_t i = 0; i < n / 2; i++)
        expTable[i] = std::polar(1.0, (inverse ? 2 : -2) * M_PI * i / n);

    // Bit-reversed addressing permutation
    for (size_t i = 0; i < n; i++) {
        size_t j = reverseBits(i, levels);
        if (j > i)
            std::swap(vec[i], vec[j]);
    }

    // Cooley-Tukey decimation-in-time radix-2 FFT
    for (size_t size = 2; size <= n; size *= 2) {
        size_t halfsize = size / 2;
        size_t tablestep = n / size;
        for (size_t i = 0; i < n; i += size) {
            for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
                Complex temp = vec[j + halfsize] * expTable[k];
                vec[j + halfsize] = vec[j] - temp;
                vec[j] += temp;
            }
        }
        if (size == n)  // Prevent overflow in 'size *= 2'
            break;
    }
}

size_t WAVFile::reverseBits(size_t val, int width) {
    size_t result = 0;
    for (int i = 0; i < width; i++, val >>= 1)
        result = (result << 1) | (val & 1U);
    return result;
}

void WAVFile::createGraph()
{
    std::ofstream outLeft("leftdata.txt");
    std::ofstream outRight("rightData.txt");

    for(size_t i = 0; i < 20000; i++)
    {
        outLeft << i / 10 << " " << leftData[i].real() << " " << leftData[i].imag() << std::endl;
        outRight << i / 10 << " " << rightData[i].real() << " " << rightData[i].imag() << std::endl;

    }
    outLeft.close();
    outRight.close();

    std::system("gnuplot -persist -e \"plot 'leftData.txt' using 1:2 title 'Real' with lines, 'leftData.txt' using 1:3 title 'Imag' with lines\"");
    std::system("gnuplot -persist -e \"plot 'rightData.txt' using 1:2 title 'Real' with lines, 'rightData.txt' using 1:3 title 'Imag' with lines\"");
}
