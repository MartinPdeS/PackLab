#include "discrete.h"


// ================================================================
// Discrete
// ================================================================

Discrete::Discrete(std::vector<double> _radii, std::vector<double> weights)
: radii(std::move(_radii))
{
    number_of_bins_ = radii.size();
    if (radii.empty())
        throw std::invalid_argument("radii must not be empty.");

    if (weights.size() != radii.size())
        throw std::invalid_argument("weights must match radii size.");

    for (double r : radii) {
        if (r <= 0.0)
            throw std::invalid_argument("radii must be positive.");

        maximum_radius = std::max(maximum_radius, r);
    }

    double weight_sum = 0.0;
    for (double w : weights) {
        if (w < 0.0)
            throw std::invalid_argument("weights must be >= 0.");
        weight_sum += w;
    }

    if (weight_sum <= 0.0)
        throw std::invalid_argument("weights must sum to a positive value.");

    // NEW: store normalized weights for to_bins()
    weights_ = std::move(weights);
    for (double& w : weights_) {
        w /= weight_sum;
    }

    // Build CDF from normalized weights_
    cumulative_probability.resize(weights_.size());
    double cumulative = 0.0;

    for (std::size_t i = 0; i < weights_.size(); ++i) {
        cumulative += weights_[i];
        cumulative_probability[i] = cumulative;
    }

    cumulative_probability.back() = 1.0;

    this->set_number_of_bins(this->number_of_bins_);

    if (number_of_bins_ > 0) {
        bin_edges_.resize(number_of_bins_ + 1);

        const double min_r = *std::min_element(radii.begin(), radii.end());
        const double max_r = maximum_radius;

        const double dr = (max_r - min_r) / static_cast<double>(number_of_bins_);

        for (std::size_t i = 0; i <= number_of_bins_; ++i)
            bin_edges_[i] = min_r + dr * static_cast<double>(i);
    }
    this->validate_bin_edges();
}


double Discrete::sample_radius(std::mt19937_64& random_generator)
{
    std::uniform_real_distribution<double> uniform_01(0.0, 1.0);
    const double u = uniform_01(random_generator);

    auto it = std::lower_bound(cumulative_probability.begin(), cumulative_probability.end(), u);

    std::size_t index = std::distance(cumulative_probability.begin(), it);
    if (index >= radii.size())
        index = radii.size() - 1;

    return radii[index];
}

int Discrete::bin_index(double r) const
{
    if (number_of_bins_ == 0)
        return -1;

    auto it = std::upper_bound(bin_edges_.begin(), bin_edges_.end(), r);
    int idx = static_cast<int>(std::distance(bin_edges_.begin(), it)) - 1;

    if (idx < 0)
        idx = 0;
    if (idx >= static_cast<int>(number_of_bins_))
        idx = static_cast<int>(number_of_bins_ - 1);

    return idx;
}


std::pair<std::vector<double>, std::vector<double>>
Discrete::to_bins() const
{
    if (radii.empty()) {
        throw std::runtime_error("Discrete::to_bins: radii is empty.");
    }
    if (weights_.size() != radii.size()) {
        throw std::runtime_error("Discrete::to_bins: weights and radii size mismatch.");
    }
    return {radii, weights_};
}
