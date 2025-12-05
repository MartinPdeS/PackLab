#pragma once

#include <vector>   // for std::vector
#include <cmath>    // for std::exp
#include <algorithm> // for std::min
#include <cctype>
#include <stdexcept>
#include <string>

#include "../utils/utils.h"
#include "../domain/domain.h"
#include "../rsa/spatial_grid_index.h"

class Result {
public:
    Domain domain;

    Result() = default;

    /*
    Constructor to initialize the Result with particle positions, radii, and domain box.
    @param particle_positions Vector of particle center positions.
    @param particle_radii Vector of particle radii.
    @param domain Simulation domain box.
    */
    Result(std::vector<Vector3d> particle_positions, std::vector<double> particle_radii, Domain domain)
        : domain(std::move(domain)), particle_positions_(std::move(particle_positions)), particle_radii_(std::move(particle_radii))
    {}

    /*
    Constructor to initialize the Result with raw particle positions, radii, and domain box.
    @param raw_positions Vector of vectors representing particle center positions.
    @param particle_radii Vector of particle radii.
    @param domain Simulation domain box.
    */
    Result(const std::vector<std::vector<double>>& raw_positions, std::vector<double> particle_radii, Domain _domain)
        : domain(std::move(_domain)), particle_radii_(std::move(particle_radii))
    {
        particle_positions_.reserve(raw_positions.size());

        for (const auto& row : raw_positions) {
            if (row.size() != 3) {
                throw std::invalid_argument(
                    "Each position entry must contain exactly three values."
                );
            }
            particle_positions_.emplace_back(row[0], row[1], row[2]);
        }
    }

    /*
    Accessors of particle positions and radii.
    @return References to the internal vectors of positions and radii.
    */
    const std::vector<Vector3d>& particle_positions() const {
        return particle_positions_;
    }

    /*
    Accessor for particle radii.
    @return Reference to the internal vector of particle radii.
    */
    const std::vector<double>& particle_radii() const {
        return particle_radii_;
    }

    /*
    Accessors for pair correlation function data.
    @return References to the internal vectors of bin centers and g(r) values.
    */
    const std::vector<double>& pair_correlation_centers() const {
        return pair_correlation_centers_;
    }

    /*
    Accessor for pair correlation values.
    @return Reference to the internal vector of g(r) values.
    */
    const std::vector<double>& pair_correlation_values() const {
        return pair_correlation_values_;
    }

    /*
    Accessor for pair correlation mean values.
    @return Reference to the internal vector of mean g(r) values.
    */
    const std::vector<double>& pair_correlation_mean_values() const {
        return pair_correlation_mean_values_;
    }

    /*
    Accessor for pair correlation standard deviation values.
    @return Reference to the internal vector of g(r) standard deviation values.
    */
    const std::vector<double>& pair_correlation_std_values() const {
        return pair_correlation_std_values_;
    }

    /*
    Compute the pair correlation function (g(r)) of the particle configuration.
    @param bins Number of bins to use for the histogram.
    @param maximum_pairs Maximum number of random pairs to sample for the calculation.
    @param method Method to use for computation: "grid" or "monte-carlo".
    @param random_seed Seed for the random number generator used to sample pairs.
    @param maximum_distance Maximum distance to consider for g(r). If zero or negative, it defaults to half the smallest box length.
    @param grid_cell_size Size of the grid cells used for spatial partitioning in the calculation.
    */
    void compute_pair_correlation_function(
        std::size_t bins,
        double maximum_distance,
        const std::string& method,
        std::size_t maximum_number_of_pairs,
        std::uint64_t random_seed,
        double grid_cell_size
    );



private:

    /*
    Fill the pair correlation histogram using Monte Carlo sampling of particle pairs.
    @param histogram Reference to the histogram vector to fill.
    @param bins Number of bins in the histogram.
    @param maximum_number_of_pairs Maximum number of random pairs to sample.
    @param random_seed Seed for the random number generator.
    @param maximum_distance Maximum distance to consider for g(r).
    @param radial_bin_width Width of each radial bin.
    */
    void fill_pair_correlation_histogram_monte_carlo(
        std::vector<std::size_t>& histogram,
        std::size_t bins,
        std::size_t maximum_number_of_pairs,
        std::uint64_t random_seed,
        double maximum_distance,
        double radial_bin_width) const;

    /*
    Compute pair correlation values from the histogram.
    @param bins Number of bins in the histogram.
    @param maximum_number_of_pairs Maximum number of random pairs sampled.
    @param random_seed Seed for the random number generator.
    @param maximum_distance Maximum distance considered for g(r).
    @return Vector of computed g(r) values.
    */
    std::vector<double> compute_pair_correlation_values_once(
        std::size_t bins,
        std::size_t maximum_number_of_pairs,
        std::uint64_t random_seed,
        double maximum_distance
    ) const;


    /*
    Fill the pair correlation histogram using a spatial grid index for efficiency.
    @param histogram Reference to the histogram vector to fill.
    @param bins Number of bins in the histogram.
    @param maximum_distance Maximum distance to consider for g(r).
    @param radial_bin_width Width of each radial bin.
    @param grid_cell_size Size of the grid cells used for spatial partitioning.
    */
    void fill_pair_correlation_histogram_grid(
        std::vector<std::size_t>& histogram,
        std::size_t bins,
        double maximum_distance,
        double radial_bin_width,
        double grid_cell_size
    ) const;

    /*
    Compute pair correlation values from the histogram.
    @param pair_correlation_values Reference to the vector to store g(r) values.
    @param histogram Reference to the histogram vector.
    @param number_of_particles Number of particles in the configuration.
    @param bins Number of bins in the histogram.
    @param radial_bin_width Width of each radial bin.
    @param maximum_distance Maximum distance considered for g(r).
    @param domain_volume Volume of the simulation domain.
    @param total_pairs_used Total number of particle pairs sampled.
    */
    std::vector<double> compute_pair_correlation_values_once(
        std::size_t bins,
        std::size_t maximum_number_of_pairs,
        std::uint64_t random_seed,
        double distance_maximum
    );

    /*
    Compute pair correlation values using a spatial grid index for efficiency.
    @param bins Number of bins in the histogram.
    @param maximum_distance Maximum distance to consider for g(r).
    @param grid_cell_size Size of the grid cells used for spatial partitioning.
    @return Vector of computed g(r) values.
    */
    std::vector<double> compute_pair_correlation_values_once_grid(
        std::size_t bins,
        double maximum_distance,
        double grid_cell_size
    ) const;

    Vector3d apply_minimum_image(Vector3d d, const Domain& domain) const;
    double sample_distance(std::size_t number_of_particles, const std::vector<Vector3d>& positions, std::mt19937_64& rng, const Domain& domain) const;

    std::vector<Vector3d> particle_positions_;
    std::vector<double> particle_radii_;

    std::vector<double> pair_correlation_centers_;
    std::vector<double> pair_correlation_values_;

    std::vector<double> pair_correlation_mean_values_;
    std::vector<double> pair_correlation_std_values_;
};
