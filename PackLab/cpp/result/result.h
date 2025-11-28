#pragma once

#include <vector>   // for std::vector
#include <cmath>    // for std::exp
#include <algorithm> // for std::min

#include "../utils/utils.h"
#include "../domain/domain.h"

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
    Result(const std::vector<std::vector<double>>& raw_positions, std::vector<double> particle_radii, Domain domain)
        : domain(std::move(domain)), particle_radii_(std::move(particle_radii))
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
    Compute the pair correlation function (g(r)) of the particle configuration.
    @param bins Number of bins to use for the histogram.
    @param maximum_pairs Maximum number of random pairs to sample for the calculation.
    @param random_seed Seed for the random number generator used to sample pairs.
    */
    void compute_pair_correlation_function(std::size_t bins, std::size_t maximum_pairs, std::uint64_t random_seed);


private:
    Vector3d apply_minimum_image(Vector3d d, const Domain& domain) const;
    double sample_distance(std::size_t number_of_particles, const std::vector<Vector3d>& positions, std::mt19937_64& rng, const Domain& domain);

    std::vector<Vector3d> particle_positions_;
    std::vector<double> particle_radii_;

    std::vector<double> pair_correlation_centers_;
    std::vector<double> pair_correlation_values_;
};
