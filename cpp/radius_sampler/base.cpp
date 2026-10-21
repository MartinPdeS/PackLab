#include "./base.h"


std::vector<double>
RadiusSampler::normalize_weights(std::vector<double> weights) const {
    double sum = 0.0;
    for (double& w : weights) {
        if (!std::isfinite(w) || w < 0.0) {
            w = 0.0;
        }
        sum += w;
    }
    if (!(sum > 0.0)) {
        throw std::runtime_error("Distribution produced zero total weight.");
    }
    for (double& w : weights) {
        w /= sum;
    }
    return weights;
}

std::pair<std::vector<double>, std::vector<double>>
RadiusSampler::edges_to_centers_and_widths_linear(const std::vector<double>& edges) const {
    if (edges.size() < 2) {
        throw std::runtime_error("bin_edges must contain at least two points.");
    }
    const std::size_t bins = edges.size() - 1;
    std::vector<double> centers(bins);
    std::vector<double> widths(bins);

    for (std::size_t i = 0; i < bins; ++i) {
        centers[i] = 0.5 * (edges[i] + edges[i + 1]);
        widths[i]  = edges[i + 1] - edges[i];
    }
    return {centers, widths};
}

















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
