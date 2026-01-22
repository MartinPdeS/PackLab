// rsa.hpp
#pragma once

#include "base.h"

class Discrete final : public RadiusSampler {
public:
    Discrete(std::vector<double> radii, std::vector<double> weights);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return maximum_radius; }
    int bin_index(double r) const override;

    std::pair<std::vector<double>, std::vector<double>> to_bins() const override;

private:
    std::vector<double> radii;
    std::vector<double> cumulative_probability;
    std::vector<double> weights_; // store raw normalized weights for to_bins
    double maximum_radius = 0.0;
};
