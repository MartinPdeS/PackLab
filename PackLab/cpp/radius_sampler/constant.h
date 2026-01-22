// rsa.hpp
#pragma once

#include "base.h"
class Constant final : public RadiusSampler {
public:
    explicit Constant(double radius, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return radius_value_; }
    int bin_index(double r) const override;

    std::pair<std::vector<double>, std::vector<double>> to_bins() const override;

private:
    double radius_value_ = 1.0;
};
