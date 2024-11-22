//
// Created by liamg on 11/22/2024.
//
#include <cmath>
#include <complex>

template <typename T, typename K>
T fxn(const T time, const K freq)
{
    return getIntensity(time) * pow(exp(1),-2.0 * M_PI * std::complex<double>(0,1) * freq * time);
}

long double trapezoidalRule(long double a, long double b, int n, double freq)
{
    long double h = (b - a) / n;
    long double sum = 0.5 * (fxn(a, freq) + fxn(b, freq));

    for (int i = 1; i < n; i++)
    {
        double x = a + i * h;
        sum += fxn(x, freq);
    }
    return h * sum;
}