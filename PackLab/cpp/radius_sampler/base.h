// rsa.hpp
#pragma once

#include <random>       // std::mt19937_64
#include <vector>       // std::vector
#include <functional>   // std::function
#include <stdexcept>
#include <utility>

#include "monte_carlo/utils/utils.h"


class RadiusSampler {
public:
    virtual ~RadiusSampler() = default;

    virtual double sample_radius(std::mt19937_64& random_generator) = 0;
    virtual double maximum_possible_radius() const = 0;
    virtual int bin_index(double radius) const = 0;

    std::size_t number_of_bins() const { return number_of_bins_; }
    void set_number_of_bins(std::size_t bins);

    virtual std::pair<std::vector<double>, std::vector<double>> to_bins() const = 0;

protected:
    double apply_binning(double value) const;

protected:
    std::size_t number_of_bins_ = 0;
    mutable std::vector<double> bin_edges_;
    void validate_bin_edges() const;

public:
    std::vector<double> normalize_weights(std::vector<double> weights) const;

    std::pair<std::vector<double>, std::vector<double>>
    edges_to_centers_and_widths_linear(const std::vector<double>& edges) const;

};







