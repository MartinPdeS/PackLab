#include "./normal.h"


Normal::Normal(double mean, double sigma, double maximum_radius_clip, int bins)
: mean_value_(mean),
  sigma_value_(sigma),
  maximum_radius_clip_value_(maximum_radius_clip)
{
    if (!std::isfinite(mean_value_) || !std::isfinite(sigma_value_) || !std::isfinite(maximum_radius_clip_value_)) {
        throw std::invalid_argument("Normal: mean, sigma, and maximum_radius_clip must be finite.");
    }

    if (mean_value_ < 0.0 || (sigma_value_ == 0.0 && mean_value_ == 0.0)) {
        throw std::invalid_argument("Normal: mean must be non-negative, and positive when sigma is zero.");
    }

    if (sigma_value_ < 0.0) {
        throw std::invalid_argument("Normal: sigma must be >= 0.");
    }

    if (maximum_radius_clip_value_ <= 0.0) {
        throw std::invalid_argument("Normal: maximum_radius_clip must be positive.");
    }

    set_number_of_bins(bins);

    if (number_of_bins_ > 0) {
        const double min_r = 0.0;
        const double max_r = maximum_radius_clip_value_;

        bin_edges_.resize(number_of_bins_ + 1);
        const double dr = (max_r - min_r) / static_cast<double>(number_of_bins_);

        for (std::size_t i = 0; i <= number_of_bins_; ++i) {
            bin_edges_[i] = min_r + dr * static_cast<double>(i);
        }
    }

    this->validate_bin_edges();
}

double Normal::sample_radius(std::mt19937_64& random_generator)
{
    if (sigma_value_ == 0.0) {
        double r = mean_value_;
        if (r > maximum_radius_clip_value_) r = maximum_radius_clip_value_;
        return apply_binning(r);
    }

    std::normal_distribution<double> dist(mean_value_, sigma_value_);

    double r;
    do {
        r = dist(random_generator);
    } while (r <= 0.0);

    if (r > maximum_radius_clip_value_) {
        r = maximum_radius_clip_value_;
    }

    return apply_binning(r);
}

int Normal::bin_index(double r) const
{
    if (number_of_bins_ == 0) {
        return 0;
    }

    auto it = std::upper_bound(bin_edges_.begin(), bin_edges_.end(), r);
    int idx = static_cast<int>(std::distance(bin_edges_.begin(), it)) - 1;

    if (idx < 0) idx = 0;
    if (idx >= static_cast<int>(number_of_bins_)) {
        idx = static_cast<int>(number_of_bins_ - 1);
    }

    return idx;
}

std::pair<std::vector<double>, std::vector<double>>
Normal::to_bins() const
{
    if (number_of_bins_ == 0) {
        throw std::runtime_error("Normal::to_bins requires bins > 0.");
    }

    validate_bin_edges();

    auto [centers, widths] = edges_to_centers_and_widths_linear(bin_edges_);
    std::vector<double> weights(centers.size(), 0.0);

    const double mu = mean_value_;
    const double sigma = sigma_value_;
    const double maximum_radius = maximum_radius_clip_value_;

    if (sigma == 0.0) {
        double radius = mu;
        if (radius > maximum_radius) {
            radius = maximum_radius;
        }

        const int index = bin_index(radius);
        if (index >= 0) {
            weights[static_cast<std::size_t>(index)] = 1.0;
        }

        return {std::move(centers), std::move(weights)};
    }

    const auto survival = [mu, sigma](double x) {
        return 0.5 * std::erfc((x - mu) / (sigma * std::sqrt(2.0)));
    };
    const double positive_mass = survival(0.0);
    for (std::size_t i = 0; i < number_of_bins_; ++i) {
        weights[i] = (survival(bin_edges_[i]) - survival(bin_edges_[i + 1])) / positive_mass;
    }
    // The positive conditional tail is clipped to the upper endpoint.
    weights.back() += survival(maximum_radius) / positive_mass;

    weights = normalize_weights(std::move(weights));
    return {std::move(centers), std::move(weights)};
}
