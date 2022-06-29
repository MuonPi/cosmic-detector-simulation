/*
SampledDistribution class for defining statistical distributions for any given
cumulative distribution function (CDF) to a desired probability distribution function.
Code modified from https://stackoverflow.com/a/58633161, original author: user DarioP
*/
#pragma once

#include <algorithm>
#include <random>
#include <stdexcept>
#include <vector>

template <typename T = double, bool Interpolate = true>
class SampledDistribution {
    struct Sample {
        T prob, value;
        Sample(const T p, const T v)
            : prob(p)
            , value(v)
        {
        }
        friend bool operator<(T p, const Sample& s) { return p < s.prob; }
    };

    std::vector<Sample> SampledCDF;

public:
    struct InvalidBounds : std::runtime_error {
        using std::runtime_error::runtime_error;
    };
    struct CDFNotMonotonic : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    template <typename F>
    SampledDistribution(F&& cdfFunc, const T low, const T high, const unsigned resolution = 1024)
    {
        if (low >= high)
            throw InvalidBounds("");
        SampledCDF.reserve(resolution);
        const T cdfLow = cdfFunc(low);
        const T cdfHigh = cdfFunc(high);
        for (unsigned i = 0; i < resolution; ++i) {
            const T x = (high - low) * i / (resolution - 1) + low;
            const T p = (cdfFunc(x) - cdfLow) / (cdfHigh - cdfLow); // normalising
            //if (p < SampledCDF.back()) throw CDFNotMonotonic("");
            SampledCDF.emplace_back(p, x);
        }
    }

    template <typename Engine>
    T operator()(Engine& g)
    {
        const T cdf = std::uniform_real_distribution<T> { 0., 1. }(g);
        auto s = std::upper_bound(SampledCDF.begin(), SampledCDF.end(), cdf);
        if (Interpolate && s != SampledCDF.begin()) {
            auto bs = s - 1;
            const T r = (cdf - bs->prob) / (s->prob - bs->prob);
            return r * bs->value + (1 - r) * s->value;
        }
        return s->value;
    }
};
