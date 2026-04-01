#include "./normal.h"


Normal::Normal(double mean, double sigma, double maximum_radius_clip, int bins)
: mean_value_(mean),
  sigma_value_(sigma),
  maximum_radius_clip_value_(maximum_radius_clip)
{
    if (!std::isfinite(mean_value_) || !std::isfinite(sigma_value_) || !std::isfinite(maximum_radius_clip_value_)) {
        throw std::invalid_argument("Normal: mean, sigma, and maximum_radius_clip must be finite.");
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
        if (r < 0.0) r = 0.0;
        if (r > maximum_radius_clip_value_) r = maximum_radius_clip_value_;
        return apply_binning(r);
    }

    std::normal_distribution<double> dist(mean_value_, sigma_value_);

    double r = dist(random_generator);

    if (r < 0.0) {
        r = 0.0;
    }

    if (r > maximum_radius_clip_value_) {
        r = maximum_radius_clip_value_;
    }

    return apply_binning(r);
}

int Normal::bin_index(double r) const
{
    if (number_of_bins_ == 0) {
        return -1;
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
        if (radius < 0.0) {
            radius = 0.0;
        }
        if (radius > maximum_radius) {
            radius = maximum_radius;
        }

        const int index = bin_index(radius);
        if (index >= 0) {
            weights[static_cast<std::size_t>(index)] = 1.0;
        }

        return {std::move(centers), std::move(weights)};
    }

    const auto normal_cdf = [mu, sigma](double x) -> double {
        return 0.5 * (1.0 + std::erf((x - mu) / (sigma * std::sqrt(2.0))));
    };

    for (std::size_t i = 0; i < number_of_bins_; ++i) {
        const double left_edge = bin_edges_[i];
        const double right_edge = bin_edges_[i + 1];

        const double clipped_left_edge = std::max(0.0, left_edge);
        const double clipped_right_edge = std::min(maximum_radius, right_edge);

        double probability_mass = 0.0;

        if (clipped_right_edge > clipped_left_edge) {
            probability_mass =
                normal_cdf(clipped_right_edge) - normal_cdf(clipped_left_edge);
        }

        weights[i] = probability_mass;
    }

    // Add clipped mass from x < 0 to the first bin
    weights.front() += normal_cdf(0.0);

    // Add clipped mass from x > maximum_radius to the last bin
    weights.back() += 1.0 - normal_cdf(maximum_radius);

    weights = normalize_weights(std::move(weights));
    return {std::move(centers), std::move(weights)};
}
