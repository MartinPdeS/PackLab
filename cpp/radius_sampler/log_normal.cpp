#include "./log_normal.h"
#include <algorithm>
#include <random>
#include <cmath>


// ================================================================
// LogNormal
// ================================================================

LogNormal::LogNormal(
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

double LogNormal::sample_radius(std::mt19937_64& random_generator)
{
    std::normal_distribution<double> stdnorm(0.0, 1.0);

    double r = std::exp(mu_value_ + sigma_value_ * stdnorm(random_generator));

    if (r > maximum_radius_clip_value_)
        r = maximum_radius_clip_value_;

    return apply_binning(r);
}


int LogNormal::bin_index(double r) const
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

std::pair<std::vector<double>, std::vector<double>>
LogNormal::to_bins() const
{
    if (number_of_bins_ == 0) {
        throw std::runtime_error("LogNormal::to_bins requires bins > 0.");
    }

    if (sigma_value_ <= 0.0) {
        throw std::runtime_error("LogNormal::to_bins requires sigma > 0.");
    }

    validate_bin_edges();

    auto [centers, widths] = edges_to_centers_and_widths_linear(bin_edges_);

    std::vector<double> weights(centers.size(), 0.0);

    const double mu = mu_value_;
    const double sigma = sigma_value_;
    const double normalizer = sigma * std::sqrt(2.0 * 3.14159265358979323846);

    for (std::size_t i = 0; i < centers.size(); ++i) {
        const double x = centers[i];
        if (x <= 0.0) {
            weights[i] = 0.0;
            continue;
        }
        const double z = (std::log(x) - mu) / sigma;
        const double pdf = std::exp(-0.5 * z * z) / (x * normalizer);
        weights[i] = pdf * widths[i];
    }

    weights = normalize_weights(std::move(weights));
    return {std::move(centers), std::move(weights)};
}

