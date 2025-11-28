// rsa.hpp
#pragma once

#include <random> // for std::mt19937_64
#include <vector> // for std::vector
#include <functional> // for std::function
#include <stdexcept>

#include "../utils/utils.h"


class RadiusSampler {
    public:
        virtual ~RadiusSampler() = default;
        virtual double sample_radius(std::mt19937_64& random_generator) = 0;
        virtual double maximum_possible_radius() const = 0;
};

class ConstantRadiusSampler final : public RadiusSampler {
    public:
        explicit ConstantRadiusSampler(double radius);

        double sample_radius(std::mt19937_64& random_generator) override;
        double maximum_possible_radius() const override;

    private:
        double radius_value_ = 1.0;
};

class UniformRadiusSampler final : public RadiusSampler {
    public:
        UniformRadiusSampler(double minimum_radius, double maximum_radius);

        double sample_radius(std::mt19937_64& random_generator) override;
        double maximum_possible_radius() const override;

    private:
        double minimum_radius_value_ = 1.0;
        double maximum_radius_value_ = 1.0;
};

class LogNormalRadiusSampler final : public RadiusSampler {
    public:
        // radius = exp(mu + sigma * N(0,1))
        LogNormalRadiusSampler(double mu, double sigma, double maximum_radius_clip);

        double sample_radius(std::mt19937_64& random_generator) override;
        double maximum_possible_radius() const override;

    private:
        double mu_value_ = 0.0;
        double sigma_value_ = 0.0;
        double maximum_radius_clip_value_ = 1.0;
};

class DiscreteRadiusSampler final : public RadiusSampler {
    public:
        // radii and weights must have same size, weights do not need to be normalized
        DiscreteRadiusSampler(std::vector<double> radii, std::vector<double> weights);

        double sample_radius(std::mt19937_64& random_generator) override;
        double maximum_possible_radius() const override;

    private:
        std::vector<double> radii_values_;
        std::vector<double> cumulative_probability_;
        double maximum_radius_value_ = 0.0;
};

class PythonCallableRadiusSampler final : public RadiusSampler {
    public:
        explicit PythonCallableRadiusSampler(std::function<double()> python_callable, double maximum_possible_radius);

        double sample_radius(std::mt19937_64& random_generator) override;
        double maximum_possible_radius() const override;

    private:
        std::function<double()> python_callable_;
        double maximum_possible_radius_value_ = 1.0;
};
