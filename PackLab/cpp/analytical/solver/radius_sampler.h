// radius_sampler.h
#pragma once

#include <cstddef>     // std::size_t
#include <random>      // std::mt19937_64
#include <string>      // std::string
#include <utility>     // std::pair
#include <vector>      // std::vector

#include "monte_carlo/utils/utils.h"

// Base class: RadiusSampler
class RadiusSampler {
public:
    virtual ~RadiusSampler() = default;

    // Core API
    virtual double sample_radius(std::mt19937_64& random_generator) = 0;
    virtual double maximum_possible_radius() const = 0;

    // Binning API for RSA sampling quantization
    virtual int bin_index(double radius) const = 0;

    std::size_t number_of_bins() const { return number_of_bins_; }
    void set_number_of_bins(std::size_t bins);

    // Discretization API for mixture solvers (meters only)
    // Returns: (bin_centers_meters, normalized_weights)
    virtual std::pair<std::vector<double>, std::vector<double>> to_bins() const = 0;

    // Discretization spacing selection
    // Valid values: "linear", "log"
    void set_bin_spacing(std::string bin_spacing);
    const std::string& bin_spacing() const { return bin_spacing_; }

protected:
    double apply_binning(double value) const;
    void validate_bin_edges() const;

protected:
    std::size_t number_of_bins_ = 0;
    mutable std::vector<double> bin_edges_; // size = number_of_bins_ + 1 when enabled
    std::string bin_spacing_ = "linear";    // used by to_bins
};

// ConstantRadiusSampler
class ConstantRadiusSampler final : public RadiusSampler {
public:
    explicit ConstantRadiusSampler(double radius, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return radius_value_; }
    int bin_index(double r) const override;

    std::pair<std::vector<double>, std::vector<double>> to_bins() const override;

private:
    double radius_value_ = 1.0;
};

// UniformRadiusSampler
class UniformRadiusSampler final : public RadiusSampler {
public:
    UniformRadiusSampler(double minimum_radius, double maximum_radius, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return maximum_radius_; }
    int bin_index(double r) const override;

    std::pair<std::vector<double>, std::vector<double>> to_bins() const override;

private:
    double minimum_radius_ = 0.0;
    double maximum_radius_ = 0.0;
};

// LogNormalRadiusSampler
class LogNormalRadiusSampler final : public RadiusSampler {
public:
    LogNormalRadiusSampler(double mu, double sigma, double maximum_radius_clip, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return maximum_radius_clip_value_; }
    int bin_index(double r) const override;

    std::pair<std::vector<double>, std::vector<double>> to_bins() const override;

private:
    double mu_value_ = 0.0;                     // ln radius mean
    double sigma_value_ = 0.0;                  // ln radius std
    double maximum_radius_clip_value_ = 1.0;    // applied in sampling
};

// DiscreteRadiusSampler
class DiscreteRadiusSampler final : public RadiusSampler {
public:
    DiscreteRadiusSampler(std::vector<double> radii, std::vector<double> weights, int bins = 0);

    double sample_radius(std::mt19937_64& random_generator) override;
    double maximum_possible_radius() const override { return maximum_radius_; }
    int bin_index(double r) const override;

    std::pair<std::vector<double>, std::vector<double>> to_bins() const override;

private:
    std::vector<double> radii_;                       // sorted ascending
    std::vector<double> weights_;                     // normalized, same order as radii_
    std::vector<double> cumulative_probability_;      // normalized CDF, same order as radii_
    double maximum_radius_ = 0.0;
};
