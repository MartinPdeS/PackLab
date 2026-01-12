#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

class PolydisperseDomain {
public:
    enum class RoundingMode {
        Floor,
        Round
    };

public:
    std::vector<double> particle_radii;
    std::vector<double> number_fractions;
    double size_meters = 0.0;
    double volume_fraction = 0.0;

    PolydisperseDomain(
        double size_meters,
        std::vector<double> particle_radii,
        double volume_fraction,
        std::vector<double> number_fractions,
        RoundingMode rounding_mode
    );

    double volume() const;

    std::vector<double> get_particle_volumes() const;
    double get_total_particle_volume() const;

    double get_mean_particle_volume_number_weighted() const;

    std::int64_t get_number_of_particles_total() const;
    std::vector<std::int64_t> get_number_of_particles_per_radius() const;

    std::vector<double> get_particle_densities_per_radius() const;
    double get_particle_density_total() const;

    std::vector<double> get_volume_fraction_per_radius() const;

    std::vector<double> sample_particle_radii(
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

private:
    RoundingMode rounding_mode_ = RoundingMode::Floor;

public:
    void print_bins(int precision) const;
    std::string bins_table(int precision) const;

};
