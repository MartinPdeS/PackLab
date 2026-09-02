// rsa.hpp
#pragma once

#include "base.h"


class LogNormal final : public RadiusSampler {
public:
    LogNormal(double mu, double sigma, double maximum_radius_clip, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return maximum_radius_clip_value_; }
    int bin_index(double r) const override;

    std::pair<std::vector<double>, std::vector<double>> to_bins() const override;

private:
    double mu_value_ = 0.0;
    double sigma_value_ = 0.0;
    double maximum_radius_clip_value_ = 1.0;
};
