#include "./uniform.h"

Uniform::Uniform(double _minimum_radius, double _maximum_radius, int bins)
: minimum_radius(_minimum_radius), maximum_radius(_maximum_radius)
{
    if (minimum_radius <= 0.0 || maximum_radius <= 0.0)
        throw std::invalid_argument("Radii must be positive.");

    if (maximum_radius < minimum_radius)
        throw std::invalid_argument("maximum_radius must be >= minimum_radius.");

    set_number_of_bins(bins);

    if (number_of_bins_ > 0) {
        bin_edges_.resize(number_of_bins_ + 1);
        const double dr = (maximum_radius - minimum_radius) /
                          static_cast<double>(number_of_bins_);
        for (std::size_t i = 0; i <= number_of_bins_; ++i)
            bin_edges_[i] = minimum_radius + dr * static_cast<double>(i);
    }
    this->validate_bin_edges();
}

double Uniform::sample_radius(std::mt19937_64& random_generator)
{
    std::uniform_real_distribution<double> dist(minimum_radius, maximum_radius);

    const double r = dist(random_generator);
    return apply_binning(r);
}


int Uniform::bin_index(double r) const
{
    if (number_of_bins_ == 0)
        return -1;

    if (r <= minimum_radius)
        return 0;

    if (r >= maximum_radius)
        return static_cast<int>(number_of_bins_ - 1);

    const double normalized =
        (r - minimum_radius) /
        (maximum_radius - minimum_radius);

    int index = static_cast<int>(normalized * number_of_bins_);
    if (index >= static_cast<int>(number_of_bins_))
        index = static_cast<int>(number_of_bins_ - 1);

    return index;
}


std::pair<std::vector<double>, std::vector<double>>
Uniform::to_bins() const
{
    if (number_of_bins_ == 0) {
        throw std::runtime_error("Uniform::to_bins requires bins > 0.");
    }

    validate_bin_edges();

    auto [centers, widths] = edges_to_centers_and_widths_linear(bin_edges_);
    auto weights = normalize_weights(std::move(widths));
    return {std::move(centers), std::move(weights)};
}


