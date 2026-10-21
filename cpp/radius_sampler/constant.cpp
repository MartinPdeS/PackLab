#include "./constant.h"

Constant::Constant(double radius, int bins)
: radius_value_(radius)
{
    if (radius <= 0.0)
        throw std::invalid_argument("Radius must be positive.");

    if (bins < 0)
        throw std::invalid_argument("Number of bins must be non-negative.");

    // There is only one distinct radius in a constant distribution. Treat it
    // as one class even when the Python caller leaves ``bins`` at its default
    // of zero; the class index is consumed by the simulator and result code.
    set_number_of_bins(1);

    bin_edges_.resize(2);
    bin_edges_[0] = radius;
    bin_edges_[1] = radius;

    this->validate_bin_edges();
}

double Constant::sample_radius(std::mt19937_64&)
{
    return apply_binning(radius_value_);
}

int Constant::bin_index(double) const
{
    // A constant sampler always represents exactly one radius class. Even
    // when optional radius binning is disabled, the simulator stores this
    // value as an array index for class-resolved statistics, so it must be a
    // valid non-negative class index.
    return 0;
}


std::pair<std::vector<double>, std::vector<double>>
Constant::to_bins() const
{
    return {{radius_value_}, {1.0}};
}
