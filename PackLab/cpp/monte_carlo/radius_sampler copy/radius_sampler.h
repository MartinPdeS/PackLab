// rsa.hpp
#pragma once

#include <random>       // std::mt19937_64
#include <vector>       // std::vector
#include <functional>   // std::function
#include <stdexcept>

#include "monte_carlo/utils/utils.h"


// -----------------------------------------------------------
// Base Class: RadiusSampler
// -----------------------------------------------------------
class RadiusSampler {
public:
    virtual ~RadiusSampler() = default;

    // Core API
    virtual double sample_radius(std::mt19937_64& random_generator) = 0;
    virtual double maximum_possible_radius() const = 0;

    // Binning API: derived classes decide how radii map to bins
    virtual int bin_index(double radius) const = 0;

    // Public accessor for number of bins
    std::size_t number_of_bins() const { return number_of_bins_; }

    // Enable or disable binning
    void set_number_of_bins(std::size_t bins);

protected:
    // Helper: quantize a sampled radius into the appropriate bin
    double apply_binning(double value) const;

protected:
    std::size_t number_of_bins_ = 0;
    mutable std::vector<double> bin_edges_; // size = bins+1 when enabled
    void validate_bin_edges() const;
};


// -----------------------------------------------------------
// ConstantRadiusSampler
// -----------------------------------------------------------
class ConstantRadiusSampler final : public RadiusSampler {
public:
    explicit ConstantRadiusSampler(double radius, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return radius_value_; }

    int bin_index(double r) const override;

private:
    double radius_value_ = 1.0;


};


// -----------------------------------------------------------
// UniformRadiusSampler
// -----------------------------------------------------------
class UniformRadiusSampler final : public RadiusSampler {
public:
    UniformRadiusSampler(double minimum_radius, double maximum_radius, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return maximum_radius; }

    int bin_index(double r) const override;

private:
    double minimum_radius;
    double maximum_radius;
};


// -----------------------------------------------------------
// LogNormalRadiusSampler
// -----------------------------------------------------------
class LogNormalRadiusSampler final : public RadiusSampler {
public:
    LogNormalRadiusSampler(double mu, double sigma, double maximum_radius_clip, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return maximum_radius_clip_value_; }

    int bin_index(double r) const override;

private:
    double mu_value_ = 0.0;
    double sigma_value_ = 0.0;
    double maximum_radius_clip_value_ = 1.0;
};


// -----------------------------------------------------------
// DiscreteRadiusSampler
// -----------------------------------------------------------
class DiscreteRadiusSampler final : public RadiusSampler {
public:
    DiscreteRadiusSampler(std::vector<double> radii, std::vector<double> weights);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return maximum_radius; }

    int bin_index(double r) const override;

private:
    std::vector<double> radii;
    std::vector<double> cumulative_probability;
    double maximum_radius = 0.0;
};
