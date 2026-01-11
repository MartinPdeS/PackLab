#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace radius_core {

inline constexpr double pi_value = 3.141592653589793238462643383279502884;

inline void require(bool condition, const char* message) {
    if (!condition) {
        throw std::invalid_argument(message);
    }
}

inline std::vector<double> normalize_weights(std::vector<double> weights) {
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

inline std::vector<double> linspace(double a, double b, int count) {
    require(count >= 2, "linspace requires count >= 2.");
    std::vector<double> out(static_cast<std::size_t>(count));
    const double step = (b - a) / static_cast<double>(count - 1);
    for (int i = 0; i < count; ++i) {
        out[static_cast<std::size_t>(i)] = a + step * static_cast<double>(i);
    }
    return out;
}

inline std::vector<double> logspace10(double a, double b, int count) {
    require(count >= 2, "logspace requires count >= 2.");
    std::vector<double> out(static_cast<std::size_t>(count));
    const double step = (b - a) / static_cast<double>(count - 1);
    for (int i = 0; i < count; ++i) {
        out[static_cast<std::size_t>(i)] = std::pow(10.0, a + step * static_cast<double>(i));
    }
    return out;
}

inline std::vector<double> make_bin_edges_meters(
    double radius_min_m,
    double radius_max_m,
    int number_of_bins,
    const std::string& bin_spacing
) {
    require(number_of_bins >= 1, "number_of_bins must be >= 1.");
    require(radius_max_m > radius_min_m, "radius_max must be strictly larger than radius_min.");

    if (bin_spacing == "linear") {
        return linspace(radius_min_m, radius_max_m, number_of_bins + 1);
    }

    if (bin_spacing == "log") {
        require(radius_min_m > 0.0, "radius_min must be > 0 for log spacing.");
        return logspace10(std::log10(radius_min_m), std::log10(radius_max_m), number_of_bins + 1);
    }

    throw std::invalid_argument("bin_spacing must be 'linear' or 'log'.");
}

inline std::pair<std::vector<double>, std::vector<double>> edges_to_centers_and_widths_meters(
    const std::vector<double>& edges_m,
    const std::string& bin_spacing
) {
    require(edges_m.size() >= 2, "edges must contain at least two points.");

    const std::size_t n = edges_m.size() - 1;
    std::vector<double> centers(n);
    std::vector<double> widths(n);

    if (bin_spacing == "linear") {
        for (std::size_t i = 0; i < n; ++i) {
            const double left = edges_m[i];
            const double right = edges_m[i + 1];
            centers[i] = 0.5 * (left + right);
            widths[i] = right - left;
        }
        return {centers, widths};
    }

    if (bin_spacing == "log") {
        for (std::size_t i = 0; i < n; ++i) {
            const double left = edges_m[i];
            const double right = edges_m[i + 1];
            centers[i] = std::sqrt(left * right);
            widths[i] = right - left;
        }
        return {centers, widths};
    }

    throw std::invalid_argument("bin_spacing must be 'linear' or 'log'.");
}

struct Bins {
    std::vector<double> radii_m;   // meters
    std::vector<double> weights;   // sum to 1
};

inline Bins delta(double radius_m) {
    return Bins{{radius_m}, {1.0}};
}

inline Bins uniform(double radius_min_m, double radius_max_m, int number_of_bins, const std::string& bin_spacing) {
    auto edges = make_bin_edges_meters(radius_min_m, radius_max_m, number_of_bins, bin_spacing);
    auto [centers, widths] = edges_to_centers_and_widths_meters(edges, bin_spacing);
    auto weights = normalize_weights(widths);
    return Bins{centers, weights};
}

inline Bins gaussian(
    double mean_m,
    double sigma_m,
    double radius_min_m,
    double radius_max_m,
    int number_of_bins,
    const std::string& bin_spacing
) {
    require(sigma_m > 0.0, "standard_deviation must be > 0.");

    auto edges = make_bin_edges_meters(radius_min_m, radius_max_m, number_of_bins, bin_spacing);
    auto [centers, widths] = edges_to_centers_and_widths_meters(edges, bin_spacing);

    std::vector<double> weights(centers.size());
    const double normalizer = sigma_m * std::sqrt(2.0 * pi_value);

    for (std::size_t i = 0; i < centers.size(); ++i) {
        const double x = centers[i];
        const double z = (x - mean_m) / sigma_m;
        const double pdf = std::exp(-0.5 * z * z) / normalizer;
        weights[i] = pdf * widths[i];
    }

    weights = normalize_weights(std::move(weights));
    return Bins{centers, weights};
}

inline Bins lognormal(
    double median_m,
    double geometric_standard_deviation,
    double radius_min_m,
    double radius_max_m,
    int number_of_bins,
    const std::string& bin_spacing
) {
    require(geometric_standard_deviation > 0.0, "geometric_standard_deviation must be > 0.");
    const double sigma = std::log(geometric_standard_deviation);
    require(sigma > 0.0, "geometric_standard_deviation must be > 1 to represent a non degenerate distribution.");
    require(median_m > 0.0, "Log normal requires strictly positive radii and median_radius.");

    auto edges = make_bin_edges_meters(radius_min_m, radius_max_m, number_of_bins, bin_spacing);
    auto [centers, widths] = edges_to_centers_and_widths_meters(edges, bin_spacing);

    const double mu = std::log(median_m);
    const double normalizer = sigma * std::sqrt(2.0 * pi_value);

    std::vector<double> weights(centers.size());
    for (std::size_t i = 0; i < centers.size(); ++i) {
        const double x = centers[i];
        require(x > 0.0, "Log normal requires strictly positive radii and median_radius.");
        const double z = (std::log(x) - mu) / sigma;
        const double pdf = std::exp(-0.5 * z * z) / (x * normalizer);
        weights[i] = pdf * widths[i];
    }

    weights = normalize_weights(std::move(weights));
    return Bins{centers, weights};
}

inline Bins discrete(std::vector<double> radii_m, std::vector<double> weights) {
    require(!radii_m.empty(), "Discrete distribution must contain at least one radius value.");
    require(radii_m.size() == weights.size(), "particle_radii and weights must have the same length.");
    weights = normalize_weights(std::move(weights));
    return Bins{std::move(radii_m), std::move(weights)};
}

}
