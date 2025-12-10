#include "radius_sampler.h"
#include <algorithm>
#include <random>
#include <cmath>

// ================================================================
// Base class implementation
// ================================================================
void RadiusSampler::validate_bin_edges() const
{
    if (number_of_bins_ == 0)
        return; // no binning required

    if (bin_edges_.empty())
        throw std::runtime_error(
            "RadiusSampler: binning enabled (bins > 0) but bin_edges_ is empty. "
            "Derived sampler must initialize bin_edges_ in its constructor."
        );
}


void RadiusSampler::set_number_of_bins(std::size_t bins)
{
    number_of_bins_ = bins;

    bin_edges_.clear();

    if (bins == 0)
        return; // binning disabled
}

double RadiusSampler::apply_binning(double value) const
{
    if (number_of_bins_ == 0)
        return value;

    validate_bin_edges();  // new safety check

    // Find the bin index using edges
    auto it = std::upper_bound(bin_edges_.begin(), bin_edges_.end(), value);
    std::size_t idx = std::distance(bin_edges_.begin(), it) - 1;

    if (idx >= number_of_bins_)
        idx = number_of_bins_ - 1;

    // Return midpoint of bin
    return 0.5 * (bin_edges_[idx] + bin_edges_[idx + 1]);
}


// ================================================================
// ConstantRadiusSampler
// ================================================================

ConstantRadiusSampler::ConstantRadiusSampler(double radius, int bins)
: radius_value_(radius)
{
    if (radius <= 0.0)
        throw std::invalid_argument("Radius must be positive.");

    set_number_of_bins(bins);

    if (number_of_bins_ > 0) {
        bin_edges_.resize(2);
        bin_edges_[0] = radius;
        bin_edges_[1] = radius;
    }

    this->validate_bin_edges();
}

double ConstantRadiusSampler::sample_radius(std::mt19937_64&)
{
    return apply_binning(radius_value_);
}

int ConstantRadiusSampler::bin_index(double) const
{
    if (number_of_bins_ == 0)
        return -1;

    return 0;
}



// ================================================================
// UniformRadiusSampler
// ================================================================

UniformRadiusSampler::UniformRadiusSampler(double minimum_radius,
                                           double maximum_radius,
                                           int bins)
: minimum_radius_value_(minimum_radius),
  maximum_radius_value_(maximum_radius)
{
    if (minimum_radius <= 0.0 || maximum_radius <= 0.0)
        throw std::invalid_argument("Radii must be positive.");

    if (maximum_radius < minimum_radius)
        throw std::invalid_argument("maximum_radius must be >= minimum_radius.");

    set_number_of_bins(bins);

    if (number_of_bins_ > 0) {
        bin_edges_.resize(number_of_bins_ + 1);
        const double dr = (maximum_radius_value_ - minimum_radius_value_) /
                          static_cast<double>(number_of_bins_);
        for (std::size_t i = 0; i <= number_of_bins_; ++i)
            bin_edges_[i] = minimum_radius_value_ + dr * static_cast<double>(i);
    }
    this->validate_bin_edges();
}

double UniformRadiusSampler::sample_radius(std::mt19937_64& random_generator)
{
    std::uniform_real_distribution<double> dist(minimum_radius_value_,
                                                maximum_radius_value_);

    const double r = dist(random_generator);
    return apply_binning(r);
}


int UniformRadiusSampler::bin_index(double r) const
{
    if (number_of_bins_ == 0)
        return -1;

    if (r <= minimum_radius_value_)
        return 0;

    if (r >= maximum_radius_value_)
        return static_cast<int>(number_of_bins_ - 1);

    const double normalized =
        (r - minimum_radius_value_) /
        (maximum_radius_value_ - minimum_radius_value_);

    int index = static_cast<int>(normalized * number_of_bins_);
    if (index >= static_cast<int>(number_of_bins_))
        index = static_cast<int>(number_of_bins_ - 1);

    return index;
}


// ================================================================
// LogNormalRadiusSampler
// ================================================================

LogNormalRadiusSampler::LogNormalRadiusSampler(
    double mu,
    double sigma,
    double maximum_radius_clip,
    int bins
)
: mu_value_(mu),
  sigma_value_(sigma),
  maximum_radius_clip_value_(maximum_radius_clip)
{
    if (sigma < 0.0)
        throw std::invalid_argument("sigma must be >= 0.");

    if (maximum_radius_clip <= 0.0)
        throw std::invalid_argument("maximum_radius_clip must be positive.");

    set_number_of_bins(bins);

    if (number_of_bins_ > 0) {
        const double min_r = 0.0;
        const double max_r = maximum_radius_clip_value_;

        bin_edges_.resize(number_of_bins_ + 1);
        const double dr = (max_r - min_r) / static_cast<double>(number_of_bins_);

        for (std::size_t i = 0; i <= number_of_bins_; ++i)
            bin_edges_[i] = min_r + dr * static_cast<double>(i);
    }
    this->validate_bin_edges();
}

double LogNormalRadiusSampler::sample_radius(std::mt19937_64& random_generator)
{
    std::normal_distribution<double> stdnorm(0.0, 1.0);

    double r = std::exp(mu_value_ + sigma_value_ * stdnorm(random_generator));

    if (r > maximum_radius_clip_value_)
        r = maximum_radius_clip_value_;

    return apply_binning(r);
}


int LogNormalRadiusSampler::bin_index(double r) const
{
    if (number_of_bins_ == 0)
        return -1;

    auto it = std::upper_bound(bin_edges_.begin(), bin_edges_.end(), r);
    int idx = static_cast<int>(std::distance(bin_edges_.begin(), it)) - 1;

    if (idx < 0) idx = 0;
    if (idx >= static_cast<int>(number_of_bins_))
        idx = static_cast<int>(number_of_bins_ - 1);

    return idx;
}


// ================================================================
// DiscreteRadiusSampler
// ================================================================

DiscreteRadiusSampler::DiscreteRadiusSampler(
    std::vector<double> radii,
    std::vector<double> weights,
    int bins
)
: radii_values_(std::move(radii))
{
    if (radii_values_.empty())
        throw std::invalid_argument("radii must not be empty.");

    if (weights.size() != radii_values_.size())
        throw std::invalid_argument("weights must match radii size.");

    double weight_sum = 0.0;

    for (double r : radii_values_) {
        if (r <= 0.0)
            throw std::invalid_argument("radii must be positive.");

        maximum_radius_value_ = std::max(maximum_radius_value_, r);
    }

    for (double w : weights) {
        if (w < 0.0)
            throw std::invalid_argument("weights must be >= 0.");
        weight_sum += w;
    }

    if (weight_sum <= 0.0)
        throw std::invalid_argument("weights must sum to a positive value.");

    cumulative_probability_.resize(weights.size());
    double cumulative = 0.0;

    for (std::size_t i = 0; i < weights.size(); ++i) {
        cumulative += weights[i] / weight_sum;
        cumulative_probability_[i] = cumulative;
    }

    cumulative_probability_.back() = 1.0;

    set_number_of_bins(bins);

    if (number_of_bins_ > 0) {
        bin_edges_.resize(number_of_bins_ + 1);

        const double min_r = *std::min_element(radii_values_.begin(), radii_values_.end());
        const double max_r = maximum_radius_value_;

        const double dr = (max_r - min_r) / static_cast<double>(number_of_bins_);

        for (std::size_t i = 0; i <= number_of_bins_; ++i)
            bin_edges_[i] = min_r + dr * static_cast<double>(i);
    }
    this->validate_bin_edges();
}

double DiscreteRadiusSampler::sample_radius(std::mt19937_64& random_generator)
{
    std::uniform_real_distribution<double> uniform_01(0.0, 1.0);
    const double u = uniform_01(random_generator);

    auto it = std::lower_bound(cumulative_probability_.begin(),
                               cumulative_probability_.end(), u);

    std::size_t index = std::distance(cumulative_probability_.begin(), it);
    if (index >= radii_values_.size())
        index = radii_values_.size() - 1;

    double r = radii_values_[index];
    return apply_binning(r);
}

int DiscreteRadiusSampler::bin_index(double r) const
{
    if (number_of_bins_ == 0)
        return -1;

    auto it = std::upper_bound(bin_edges_.begin(), bin_edges_.end(), r);
    int idx = static_cast<int>(std::distance(bin_edges_.begin(), it)) - 1;

    if (idx < 0) idx = 0;
    if (idx >= static_cast<int>(number_of_bins_))
        idx = static_cast<int>(number_of_bins_ - 1);

    return idx;
}
