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


class PYDomain {
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

    PYDomain(
        double size,
        std::vector<double> radii,
        double volume_fraction,
        std::vector<double> number_fractions,
        RoundingMode rounding_mode
    );

    /*
    Get the cubic volume of the domain in cubic meters.
    @returns double volume in cubic meters.
    */
    double get_volume() const;

    /*
    Get the per bin particle volumes in cubic meters.
    @returns std::vector<double> of per bin particle volumes in cubic meters.
    */
    std::vector<double> get_particle_volumes() const;

    /*
    Get the total particle volume in cubic meters.
    @returns double total particle volume in cubic meters.
    */
    double get_total_particle_volume() const;

    /*
    Get the number-weighted mean particle volume in cubic meters.
    @returns double mean particle volume in cubic meters.
    */
    double get_mean_particle_volume_number_weighted() const;

    /*
    Get the total number of particles in the domain.
    @returns std::int64_t total number of particles.
    */
    std::int64_t get_number_of_particles_total() const;

    /*
    Get the per bin number of particles in the domain.
    @returns std::vector<std::int64_t> of per bin number of particles.
    */
    std::vector<std::int64_t> get_number_of_particles_per_radius() const;

    /*
    Get the per bin particle number densities in 1/meter**3.
    @returns std::vector<double> of per bin particle number densities in 1/meter**3.
    */
    std::vector<double> get_particle_densities_per_radius() const;

    /*
    Get the total particle number density in 1/meter**3.
    @returns double total particle number density in 1/meter**3.
    */
    double get_particle_density_total() const;

    /*
    Get the per bin occupied volume fractions.
    @returns std::vector<double> of per bin occupied volume fractions (dimensionless).
    */
    std::vector<double> get_volume_fraction_per_radius() const;

    /*
    Sample particle radii according to number_fractions.
    @param number_of_samples Number of particle radii to sample.
    @param seed RNG seed.
    @returns std::vector<double> of sampled particle radii in meters.
    */
    std::vector<double> sample_radii(std::int64_t number_of_samples, std::uint64_t seed) const;

private:
    /*
    Validate and normalize number fractions.
    @param number_fractions Reference to vector of number fractions to validate and normalize.
    @param expected_size Expected size of the number_fractions vector.
    */
    static void validate_and_normalize_number_fractions(std::vector<double>& number_fractions, std::size_t expected_size);

    /*
    Apply rounding mode to a scalar value.
    @param value Value to round.
    @returns std::int64_t rounded value.
    */
    std::int64_t apply_rounding_scalar(double value) const;

    /*
    Apply rounding mode to a vector of values.
    @param values Vector of values to round.
    @returns std::vector<std::int64_t> of rounded values.
    */
    std::vector<std::int64_t> apply_rounding_vector(const std::vector<double>& values) const;

public:
    void print_bins(int precision) const;
    std::string bins_table(int precision) const;

};
