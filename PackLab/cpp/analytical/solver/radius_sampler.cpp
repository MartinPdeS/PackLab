// radius_sampler.cpp
#include "radius_sampler.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

static constexpr double pi_value = 3.141592653589793238462643383279502884;

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::invalid_argument(message);
    }
}

std::vector<double> normalize_weights(std::vector<double> weights) {
    double total = 0.0;
    for (double& w : weights) {
        if (!std::isfinite(w) || w < 0.0) {
            w = 0.0;
        }
        total += w;
    }
    require(total > 0.0, "Distribution produced zero total weight. Check parameters or range.");
    for (double& w : weights) {
        w /= total;
    }
    return weights;
}

std::vector<double> make_edges(double radius_min_m, double radius_max_m, std::size_t bins, const std::string& bin_spacing) {
    require(bins >= 1, "number_of_bins must be >= 1.");
    require(radius_max_m > radius_min_m, "radius_max must be strictly larger than radius_min.");

    std::vector<double> edges(bins + 1);

    if (bin_spacing == "linear") {
        const double step = (radius_max_m - radius_min_m) / static_cast<double>(bins);
        for (std::size_t i = 0; i <= bins; ++i) {
            edges[i] = radius_min_m + step * static_cast<double>(i);
        }
        return edges;
    }

    if (bin_spacing == "log") {
        require(radius_min_m > 0.0, "radius_min must be > 0 for log spacing.");
        const double log_min = std::log10(radius_min_m);
        const double log_max = std::log10(radius_max_m);
        const double step = (log_max - log_min) / static_cast<double>(bins);
        for (std::size_t i = 0; i <= bins; ++i) {
            edges[i] = std::pow(10.0, log_min + step * static_cast<double>(i));
        }
        return edges;
    }

    throw std::invalid_argument("bin_spacing must be 'linear' or 'log'.");
}

std::pair<std::vector<double>, std::vector<double>> edges_to_centers_widths(
    const std::vector<double>& edges,
    const std::string& bin_spacing
) {
    require(edges.size() >= 2, "edges must contain at least two points.");

    const std::size_t bins = edges.size() - 1;
    std::vector<double> centers(bins);
    std::vector<double> widths(bins);

    if (bin_spacing == "linear") {
        for (std::size_t i = 0; i < bins; ++i) {
            centers[i] = 0.5 * (edges[i] + edges[i + 1]);
            widths[i]  = edges[i + 1] - edges[i];
        }
        return {centers, widths};
    }

    if (bin_spacing == "log") {
        for (std::size_t i = 0; i < bins; ++i) {
            centers[i] = std::sqrt(edges[i] * edges[i + 1]);
            widths[i]  = edges[i + 1] - edges[i];
        }
        return {centers, widths};
    }

    throw std::invalid_argument("bin_spacing must be 'linear' or 'log'.");
}

// Acklam inverse CDF approximation
double inverse_standard_normal_cdf(double p) {
    require(p > 0.0 && p < 1.0, "inverse_standard_normal_cdf requires p in (0, 1).");

    const double a1 = -3.969683028665376e+01;
    const double a2 =  2.209460984245205e+02;
    const double a3 = -2.759285104469687e+02;
    const double a4 =  1.383577518672690e+02;
    const double a5 = -3.066479806614716e+01;
    const double a6 =  2.506628277459239e+00;

    const double b1 = -5.447609879822406e+01;
    const double b2 =  1.615858368580409e+02;
    const double b3 = -1.556989798598866e+02;
    const double b4 =  6.680131188771972e+01;
    const double b5 = -1.328068155288572e+01;

    const double c1 = -7.784894002430293e-03;
    const double c2 = -3.223964580411365e-01;
    const double c3 = -2.400758277161838e+00;
    const double c4 = -2.549732539343734e+00;
    const double c5 =  4.374664141464968e+00;
    const double c6 =  2.938163982698783e+00;

    const double d1 =  7.784695709041462e-03;
    const double d2 =  3.224671290700398e-01;
    const double d3 =  2.445134137142996e+00;
    const double d4 =  3.754408661907416e+00;

    const double p_low  = 0.02425;
    const double p_high = 1.0 - p_low;

    if (p < p_low) {
        const double q = std::sqrt(-2.0 * std::log(p));
        return (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
               ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0);
    }

    if (p > p_high) {
        const double q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
                 ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0);
    }

    const double q = p - 0.5;
    const double r = q * q;
    return (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6) * q /
           (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0);
}

} // namespace

// Base class
void RadiusSampler::set_number_of_bins(std::size_t bins) {
    number_of_bins_ = bins;
    bin_edges_.clear();
}

void RadiusSampler::set_bin_spacing(std::string bin_spacing) {
    if (bin_spacing != "linear" && bin_spacing != "log") {
        throw std::invalid_argument("bin_spacing must be 'linear' or 'log'.");
    }
    bin_spacing_ = std::move(bin_spacing);
}

void RadiusSampler::validate_bin_edges() const {
    if (number_of_bins_ == 0) {
        return;
    }

    if (bin_edges_.size() != number_of_bins_ + 1) {
        throw std::runtime_error("RadiusSampler: bin_edges_ size must be number_of_bins_ + 1.");
    }

    for (std::size_t i = 0; i + 1 < bin_edges_.size(); ++i) {
        const double left = bin_edges_[i];
        const double right = bin_edges_[i + 1];

        if (!std::isfinite(left) || !std::isfinite(right)) {
            throw std::runtime_error("RadiusSampler: bin_edges_ must be finite.");
        }
        if (right < left) {
            throw std::runtime_error("RadiusSampler: bin_edges_ must be nondecreasing.");
        }
    }
}

double RadiusSampler::apply_binning(double value) const {
    if (number_of_bins_ == 0) {
        return value;
    }

    validate_bin_edges();

    auto it = std::upper_bound(bin_edges_.begin(), bin_edges_.end(), value);

    std::size_t idx = 0;
    if (it == bin_edges_.begin()) {
        idx = 0;
    } else {
        idx = static_cast<std::size_t>(std::distance(bin_edges_.begin(), it) - 1);
        if (idx >= number_of_bins_) {
            idx = number_of_bins_ - 1;
        }
    }

    return 0.5 * (bin_edges_[idx] + bin_edges_[idx + 1]);
}

// ConstantRadiusSampler
ConstantRadiusSampler::ConstantRadiusSampler(double radius, int bins)
    : radius_value_(radius) {
    require(radius_value_ > 0.0, "Radius must be positive.");
    require(bins >= 0, "bins must be >= 0.");

    if (bins > 0) {
        // Binning is meaningless for a constant sampler, keep a single bin
        set_number_of_bins(1);
        bin_edges_.assign(2, radius_value_);
        validate_bin_edges();
    } else {
        set_number_of_bins(0);
    }
}

double ConstantRadiusSampler::sample_radius(std::mt19937_64&) {
    return apply_binning(radius_value_);
}

int ConstantRadiusSampler::bin_index(double) const {
    if (number_of_bins_ == 0) {
        return -1;
    }
    return 0;
}

std::pair<std::vector<double>, std::vector<double>> ConstantRadiusSampler::to_bins() const {
    return {{radius_value_}, {1.0}};
}

// UniformRadiusSampler
UniformRadiusSampler::UniformRadiusSampler(double minimum_radius, double maximum_radius, int bins)
    : minimum_radius_(minimum_radius), maximum_radius_(maximum_radius) {
    require(minimum_radius_ > 0.0, "Radii must be positive.");
    require(maximum_radius_ > 0.0, "Radii must be positive.");
    require(maximum_radius_ >= minimum_radius_, "maximum_radius must be >= minimum_radius.");
    require(bins >= 0, "bins must be >= 0.");

    set_number_of_bins(static_cast<std::size_t>(bins));

    // These edges are used only for RSA quantization. Use linear edges here.
    if (number_of_bins_ > 0) {
        bin_edges_.resize(number_of_bins_ + 1);
        const double dr = (maximum_radius_ - minimum_radius_) / static_cast<double>(number_of_bins_);
        for (std::size_t i = 0; i <= number_of_bins_; ++i) {
            bin_edges_[i] = minimum_radius_ + dr * static_cast<double>(i);
        }
        validate_bin_edges();
    }
}

double UniformRadiusSampler::sample_radius(std::mt19937_64& random_generator) {
    std::uniform_real_distribution<double> dist(minimum_radius_, maximum_radius_);
    return apply_binning(dist(random_generator));
}

int UniformRadiusSampler::bin_index(double r) const {
    if (number_of_bins_ == 0) {
        return -1;
    }

    if (r <= minimum_radius_) {
        return 0;
    }
    if (r >= maximum_radius_) {
        return static_cast<int>(number_of_bins_ - 1);
    }

    const double normalized = (r - minimum_radius_) / (maximum_radius_ - minimum_radius_);
    int index = static_cast<int>(normalized * static_cast<double>(number_of_bins_));
    if (index >= static_cast<int>(number_of_bins_)) {
        index = static_cast<int>(number_of_bins_ - 1);
    }
    return index;
}

std::pair<std::vector<double>, std::vector<double>> UniformRadiusSampler::to_bins() const {
    require(number_of_bins_ >= 1, "number_of_bins must be >= 1 to compute bins.");

    const std::string& spacing = bin_spacing_;
    auto edges = make_edges(minimum_radius_, maximum_radius_, number_of_bins_, spacing);
    auto centers_widths = edges_to_centers_widths(edges, spacing);

    auto weights = normalize_weights(std::move(centers_widths.second));
    return {std::move(centers_widths.first), std::move(weights)};
}

// LogNormalRadiusSampler
LogNormalRadiusSampler::LogNormalRadiusSampler(double mu, double sigma, double maximum_radius_clip, int bins)
    : mu_value_(mu),
      sigma_value_(sigma),
      maximum_radius_clip_value_(maximum_radius_clip) {
    require(sigma_value_ >= 0.0, "sigma must be >= 0.");
    require(maximum_radius_clip_value_ > 0.0, "maximum_radius_clip must be positive.");
    require(bins >= 0, "bins must be >= 0.");

    set_number_of_bins(static_cast<std::size_t>(bins));

    // These edges are used only for RSA quantization. Use linear edges from 0 to clip.
    if (number_of_bins_ > 0) {
        bin_edges_.resize(number_of_bins_ + 1);
        const double dr = maximum_radius_clip_value_ / static_cast<double>(number_of_bins_);
        for (std::size_t i = 0; i <= number_of_bins_; ++i) {
            bin_edges_[i] = dr * static_cast<double>(i);
        }
        validate_bin_edges();
    }
}

double LogNormalRadiusSampler::sample_radius(std::mt19937_64& random_generator) {
    std::normal_distribution<double> stdnorm(0.0, 1.0);

    double r = std::exp(mu_value_ + sigma_value_ * stdnorm(random_generator));

    // This sampler is clipped by construction via maximum_radius_clip_value_
    if (r > maximum_radius_clip_value_) {
        r = maximum_radius_clip_value_;
    }

    return apply_binning(r);
}

int LogNormalRadiusSampler::bin_index(double r) const {
    if (number_of_bins_ == 0) {
        return -1;
    }

    validate_bin_edges();

    auto it = std::upper_bound(bin_edges_.begin(), bin_edges_.end(), r);
    int idx = static_cast<int>(std::distance(bin_edges_.begin(), it)) - 1;

    if (idx < 0) {
        idx = 0;
    }
    if (idx >= static_cast<int>(number_of_bins_)) {
        idx = static_cast<int>(number_of_bins_ - 1);
    }
    return idx;
}

std::pair<std::vector<double>, std::vector<double>> LogNormalRadiusSampler::to_bins() const {
    require(number_of_bins_ >= 1, "number_of_bins must be >= 1 to compute bins.");
    require(sigma_value_ > 0.0, "sigma must be > 0.");

    // Numerical truncation only, no user provided min value.
    // Upper bound respects the sampler maximum_radius_clip_value_ so it matches sampling behavior.
    const double epsilon = 1e-12;
    const double p_high = 1.0 - 0.5 * epsilon;
    const double z = inverse_standard_normal_cdf(p_high);

    const double radius_min_m = std::exp(mu_value_ - z * sigma_value_);
    const double unbounded_max = std::exp(mu_value_ + z * sigma_value_);
    const double radius_max_m = std::min(unbounded_max, maximum_radius_clip_value_);

    require(radius_min_m > 0.0, "Log normal requires a strictly positive radius_min.");
    require(radius_max_m > radius_min_m, "Computed bounds are invalid. Check mu, sigma, and maximum_radius_clip.");

    const std::string& spacing = bin_spacing_;
    auto edges = make_edges(radius_min_m, radius_max_m, number_of_bins_, spacing);
    auto centers_widths = edges_to_centers_widths(edges, spacing);

    const std::vector<double>& centers = centers_widths.first;
    const std::vector<double>& widths  = centers_widths.second;

    std::vector<double> weights(centers.size());

    const double sigma = sigma_value_;
    const double mu = mu_value_;
    const double normalizer = sigma * std::sqrt(2.0 * pi_value);

    for (std::size_t i = 0; i < centers.size(); ++i) {
        const double x = centers[i];
        require(x > 0.0, "Log normal requires strictly positive radii.");
        const double z_local = (std::log(x) - mu) / sigma;
        const double pdf = std::exp(-0.5 * z_local * z_local) / (x * normalizer);
        weights[i] = pdf * widths[i];
    }

    weights = normalize_weights(std::move(weights));
    return {centers_widths.first, std::move(weights)};
}

// DiscreteRadiusSampler
DiscreteRadiusSampler::DiscreteRadiusSampler(std::vector<double> radii, std::vector<double> weights, int bins)
    : radii_(std::move(radii)) {
    require(bins >= 0, "bins must be >= 0.");
    require(!radii_.empty(), "radii must not be empty.");
    require(weights.size() == radii_.size(), "weights must match radii size.");

    for (double r : radii_) {
        require(std::isfinite(r) && r > 0.0, "radii must be finite and positive.");
        maximum_radius_ = std::max(maximum_radius_, r);
    }

    for (double w : weights) {
        require(std::isfinite(w) && w >= 0.0, "weights must be finite and >= 0.");
    }

    weights_ = normalize_weights(std::move(weights));

    // Sort by radius and permute weights accordingly
    std::vector<std::size_t> order(radii_.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return radii_[a] < radii_[b];
    });

    std::vector<double> radii_sorted(radii_.size());
    std::vector<double> weights_sorted(weights_.size());

    for (std::size_t i = 0; i < order.size(); ++i) {
        radii_sorted[i] = radii_[order[i]];
        weights_sorted[i] = weights_[order[i]];
    }

    radii_ = std::move(radii_sorted);
    weights_ = normalize_weights(std::move(weights_sorted));

    // Build CDF
    cumulative_probability_.resize(weights_.size());
    double cumulative = 0.0;
    for (std::size_t i = 0; i < weights_.size(); ++i) {
        cumulative += weights_[i];
        cumulative_probability_[i] = cumulative;
    }
    cumulative_probability_.back() = 1.0;

    // Optional RSA quantization bins
    set_number_of_bins(static_cast<std::size_t>(bins));
    if (number_of_bins_ > 0) {
        const double min_r = radii_.front();
        const double max_r = radii_.back();
        require(max_r > min_r, "Discrete sampler requires max radius > min radius when binning is enabled.");

        bin_edges_.resize(number_of_bins_ + 1);
        const double dr = (max_r - min_r) / static_cast<double>(number_of_bins_);
        for (std::size_t i = 0; i <= number_of_bins_; ++i) {
            bin_edges_[i] = min_r + dr * static_cast<double>(i);
        }
        validate_bin_edges();
    }
}

double DiscreteRadiusSampler::sample_radius(std::mt19937_64& random_generator) {
    std::uniform_real_distribution<double> uniform_01(0.0, 1.0);
    const double u = uniform_01(random_generator);

    auto it = std::lower_bound(cumulative_probability_.begin(), cumulative_probability_.end(), u);
    std::size_t index = static_cast<std::size_t>(std::distance(cumulative_probability_.begin(), it));
    if (index >= radii_.size()) {
        index = radii_.size() - 1;
    }

    return apply_binning(radii_[index]);
}

int DiscreteRadiusSampler::bin_index(double r) const {
    if (number_of_bins_ == 0) {
        return -1;
    }

    validate_bin_edges();

    auto it = std::upper_bound(bin_edges_.begin(), bin_edges_.end(), r);
    int idx = static_cast<int>(std::distance(bin_edges_.begin(), it)) - 1;

    if (idx < 0) {
        idx = 0;
    }
    if (idx >= static_cast<int>(number_of_bins_)) {
        idx = static_cast<int>(number_of_bins_ - 1);
    }
    return idx;
}

std::pair<std::vector<double>, std::vector<double>> DiscreteRadiusSampler::to_bins() const {
    require(!radii_.empty(), "DiscreteRadiusSampler::to_bins: radii must not be empty.");
    require(radii_.size() == weights_.size(), "DiscreteRadiusSampler::to_bins: radii and weights must have the same length.");
    return {radii_, weights_};
}
