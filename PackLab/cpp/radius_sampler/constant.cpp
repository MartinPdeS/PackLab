#include "./constant.h"

Constant::Constant(double radius, int bins)
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

double Constant::sample_radius(std::mt19937_64&)
{
    return apply_binning(radius_value_);
}

int Constant::bin_index(double) const
{
    if (number_of_bins_ == 0)
        return -1;

    return 0;
}


std::pair<std::vector<double>, std::vector<double>>
Constant::to_bins() const
{
    return {{radius_value_}, {1.0}};
}
