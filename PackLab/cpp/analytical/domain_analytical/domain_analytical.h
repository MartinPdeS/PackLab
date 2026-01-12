#pragma once

#include <cstddef>  // for std::size_t
#include <cstdint>  // for std::int64_t, std::uint64_t
#include <string>   // for std::string
#include <vector>   // for std::vector
#include <cmath>    // for std::round
#include <stdexcept> // for std::invalid_argument
#include <iomanip>  // for std::setprecision
#include <iostream> // for std::cout, std::endl
#include <sstream>  // for std::ostringstream
#include <random>  // for std::mt19937_64, std::random_device
#include <algorithm> // for std::sort, std::max


class Domain {
public:
    enum class RoundingMode {Floor, Round};

public:
    double size = 0.0;
    std::vector<double> radii;
    double volume_fraction = 0.0;
    std::vector<double> number_fractions;
private:
    RoundingMode rounding_mode_ = RoundingMode::Floor;
public:

    Domain(
        double size,
        std::vector<double> radii,
        double volume_fraction,
        std::vector<double> number_fractions,
        RoundingMode rounding_mode
    );

    double get_volume() const;

    std::vector<double> get_particle_volumes() const;
    double get_total_particle_volume() const;

    double get_mean_particle_volume_number_weighted() const;

    std::int64_t get_number_of_particles_total() const;
    std::vector<std::int64_t> get_number_of_particles_per_radius() const;

    std::vector<double> get_particle_densities_per_radius() const;
    double get_particle_density_total() const;

    std::vector<double> get_volume_fraction_per_radius() const;

    std::vector<double> sample_radii(
        std::int64_t number_of_samples,
        std::uint64_t seed
    ) const;

private:
    static void validate_and_normalize_number_fractions(
        std::vector<double>& number_fractions,
        std::size_t expected_size
    );

    std::int64_t apply_rounding_scalar(double value) const;
    std::vector<std::int64_t> apply_rounding_vector(const std::vector<double>& values) const;

public:
    void print_bins(int precision) const;
    std::string bins_table(int precision) const;

};
