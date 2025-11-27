// radius_sampler.cpp
#include "radius_sampler.h"

ConstantRadiusSampler::ConstantRadiusSampler(double radius)
: radius_value_(radius) {
    if (radius_value_ <= 0.0) {
        throw std::invalid_argument("Radius must be positive.");
    }
}

double ConstantRadiusSampler::sample_radius(std::mt19937_64&) { return radius_value_; }
double ConstantRadiusSampler::maximum_possible_radius() const { return radius_value_; }

UniformRadiusSampler::UniformRadiusSampler(double minimum_radius, double maximum_radius)
: minimum_radius_value_(minimum_radius)
, maximum_radius_value_(maximum_radius) {
    if (minimum_radius_value_ <= 0.0 || maximum_radius_value_ <= 0.0) {
        throw std::invalid_argument("Radii must be positive.");
    }
    if (maximum_radius_value_ < minimum_radius_value_) {
        throw std::invalid_argument("maximum_radius must be >= minimum_radius.");
    }
}

double UniformRadiusSampler::sample_radius(std::mt19937_64& random_generator) {
    std::uniform_real_distribution<double> distribution(minimum_radius_value_, maximum_radius_value_);
    return distribution(random_generator);
}

double UniformRadiusSampler::maximum_possible_radius() const { return maximum_radius_value_; }

LogNormalRadiusSampler::LogNormalRadiusSampler(double mu, double sigma, double maximum_radius_clip)
: mu_value_(mu)
, sigma_value_(sigma)
, maximum_radius_clip_value_(maximum_radius_clip) {
    if (sigma_value_ < 0.0) {
        throw std::invalid_argument("sigma must be >= 0.");
    }
    if (maximum_radius_clip_value_ <= 0.0) {
        throw std::invalid_argument("maximum_radius_clip must be positive.");
    }
}

double LogNormalRadiusSampler::sample_radius(std::mt19937_64& random_generator) {
    std::normal_distribution<double> standard_normal(0.0, 1.0);
    const double radius = std::exp(mu_value_ + sigma_value_ * standard_normal(random_generator));
    return std::min(radius, maximum_radius_clip_value_);
}

double LogNormalRadiusSampler::maximum_possible_radius() const { return maximum_radius_clip_value_; }

DiscreteRadiusSampler::DiscreteRadiusSampler(std::vector<double> radii, std::vector<double> weights)
: radii_values_(std::move(radii)) {
    if (radii_values_.empty()) {
        throw std::invalid_argument("radii must not be empty.");
    }
    if (weights.size() != radii_values_.size()) {
        throw std::invalid_argument("weights must have the same size as radii.");
    }
    for (double r : radii_values_) {
        if (r <= 0.0) throw std::invalid_argument("All radii must be positive.");
        maximum_radius_value_ = std::max(maximum_radius_value_, r);
    }

    double weight_sum = 0.0;
    for (double w : weights) {
        if (w < 0.0) throw std::invalid_argument("weights must be >= 0.");
        weight_sum += w;
    }
    if (weight_sum <= 0.0) {
        throw std::invalid_argument("weights must sum to a positive value.");
    }

    cumulative_probability_.resize(weights.size());
    double cumulative = 0.0;
    for (std::size_t index = 0; index < weights.size(); ++index) {
        cumulative += weights[index] / weight_sum;
        cumulative_probability_[index] = cumulative;
    }
    cumulative_probability_.back() = 1.0;
}

double DiscreteRadiusSampler::sample_radius(std::mt19937_64& random_generator) {
    std::uniform_real_distribution<double> uniform_01(0.0, 1.0);
    const double u = uniform_01(random_generator);

    auto iterator = std::lower_bound(cumulative_probability_.begin(), cumulative_probability_.end(), u);
    const std::size_t index = static_cast<std::size_t>(std::distance(cumulative_probability_.begin(), iterator));
    return radii_values_[std::min(index, radii_values_.size() - 1)];
}

double DiscreteRadiusSampler::maximum_possible_radius() const { return maximum_radius_value_; }

PythonCallableRadiusSampler::PythonCallableRadiusSampler(std::function<double()> python_callable, double maximum_possible_radius)
: python_callable_(std::move(python_callable))
, maximum_possible_radius_value_(maximum_possible_radius) {
    if (!python_callable_) {
        throw std::invalid_argument("python_callable must be valid.");
    }
    if (maximum_possible_radius_value_ <= 0.0) {
        throw std::invalid_argument("maximum_possible_radius must be positive.");
    }
}

double PythonCallableRadiusSampler::sample_radius(std::mt19937_64&) {
    const double value = python_callable_();
    if (!(value > 0.0)) {
        throw std::invalid_argument("PythonCallableRadiusSampler returned a non positive radius.");
    }
    return value;
}

double PythonCallableRadiusSampler::maximum_possible_radius() const { return maximum_possible_radius_value_; }
