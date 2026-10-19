// rsa.hpp
#pragma once

#include "base.h"

class Uniform final : public RadiusSampler {
public:
    Uniform(double minimum_radius, double maximum_radius, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return maximum_radius; }
    int bin_index(double r) const override;

    std::pair<std::vector<double>, std::vector<double>> to_bins() const override;

private:
    double minimum_radius;
    double maximum_radius;
};

