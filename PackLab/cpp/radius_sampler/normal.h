// rsa.hpp
#pragma once

#include "base.h"
#include <algorithm>
#include <random>
#include <cmath>


class Normal : public RadiusSampler {
public:
    Normal() = default;
    Normal(double mean, double sigma, double maximum_radius_clip, int bins);

    double sample_radius(std::mt19937_64& random_generator) override;

    int bin_index(double r) const override;

    std::pair<std::vector<double>, std::vector<double>> to_bins() const override;

    double maximum_possible_radius() const override { return 0; }

private:
    double mean_value_;
    double sigma_value_;
    double maximum_radius_clip_value_;
};
