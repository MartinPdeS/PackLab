#pragma once

#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>

#include "../utils/utils.h"
#include "../statistics/statistics.h"
#include "../domain/domain.h"
#include "../rsa/spatial_grid_index.h"
#include "../rsa/sphere_configuration.h"

/*
Result stores a completed RSA configuration along with utilities for
computing total and partial pair correlation functions g(r) and g_ij(r).
*/
class Result {
public:
    Domain domain;
    int number_of_classes_ = 0;

    /*
    Each particle is assigned a discrete class index in [0, number_of_classes_-1].
    This enables partial pair-correlation functions g_ij(r).
    */
    SphereConfiguration sphere_configuration;

    // Total g(r)
    std::vector<double> pair_correlation_centers_;
    std::vector<double> pair_correlation_values_;

    // Optional: mean and std if user samples repeatedly
    std::vector<double> pair_correlation_mean_values_;
    std::vector<double> pair_correlation_std_values_;

    Statistics statistics;

    Result() = default;


    /*
    Overload including class assignments for polydisperse systems.

    @param sphere_configuration Sphere configuration containing positions, radii, and classes.
    @param domain Simulation domain.
    @param number_of_classes Total number of distinct classes.
    */
    Result(
        SphereConfiguration sphere_configuration,
        Domain domain,
        Statistics statistics,
        int number_of_classes
    )
        : domain(std::move(domain)),
          number_of_classes_(number_of_classes),
          sphere_configuration(std::move(sphere_configuration)),
          statistics(std::move(statistics))
    {
        if (this->sphere_configuration.center_positions_values_.size() != this->sphere_configuration.radii_values_.size())
            throw std::invalid_argument("positions and radii must have same length.");

        if (!this->sphere_configuration.class_index_values_.empty() &&
            this->sphere_configuration.class_index_values_.size() != this->sphere_configuration.center_positions_values_.size())
        {
            throw std::invalid_argument(
                "classes must have same length as positions."
            );
        }
    }

    // --------------------- Accessors -------------------------

    const std::vector<Vector3d>& particle_positions() const {
        return this->sphere_configuration.center_positions_values_;
    }

    const std::vector<double>& particle_radii() const {
        return this->sphere_configuration.radii_values_;
    }

    const std::vector<double>& pair_correlation_centers() const {
        return pair_correlation_centers_;
    }

    const std::vector<double>& pair_correlation_values() const {
        return pair_correlation_values_;
    }

    const std::vector<double>& pair_correlation_mean_values() const {
        return pair_correlation_mean_values_;
    }

    const std::vector<double>& pair_correlation_std_values() const {
        return pair_correlation_std_values_;
    }

    // --------------------- Total g(r) -------------------------

    /*
    Compute total pair correlation function g(r).

    @param bins Number of distance bins.
    @param maximum_distance Maximum separation to include.
    @param method "grid" or "monte-carlo".
    @param maximum_number_of_pairs Number of sampled pairs (MC only).
    @param random_seed RNG seed.
    @param grid_cell_size Spatial grid cell size (grid method only).
    */
    void compute_pair_correlation_function(
        std::size_t bins,
        double maximum_distance,
        const std::string& method,
        std::size_t maximum_number_of_pairs,
        std::uint64_t random_seed,
        double grid_cell_size
    );

    // ---------------- Partial g_ij(r) -------------------------

    /*
    Compute partial pair correlation g_ij(r) for class-labeled particles.

    @param number_of_distance_bins Number of r-bins.
    @param maximum_pairs Maximum number of random pairs to sample.
    @param random_seed RNG seed.

    @return (bin_centers, g_ij[r][i][j]):
            - bin_centers[k] is center of bin k,
            - g_ij[k][i][j] is the partial pair-correlation between class i and j at bin k.
    */
    std::tuple<
        std::vector<double>,
        std::vector<std::vector<std::vector<double>>>
    >
    compute_partial_pair_correlation_function(
        std::size_t number_of_distance_bins,
        std::size_t maximum_pairs,
        std::uint64_t random_seed
    ) const;

private:
    // ---------------- Monte Carlo histogram fill ----------------

    void fill_pair_correlation_histogram_monte_carlo(
        std::vector<std::size_t>& histogram,
        std::size_t bins,
        std::size_t maximum_number_of_pairs,
        std::uint64_t random_seed,
        double maximum_distance,
        double radial_bin_width
    ) const;

    // ---------------- MC g(r) computation --------------------

    std::vector<double> compute_pair_correlation_values_once(
        std::size_t bins,
        std::size_t maximum_number_of_pairs,
        std::uint64_t random_seed,
        double maximum_distance
    ) const;

    // ---------------- Grid histogram fill --------------------

    void fill_pair_correlation_histogram_grid(
        std::vector<std::size_t>& histogram,
        std::size_t bins,
        double maximum_distance,
        double radial_bin_width,
        double grid_cell_size
    ) const;

    // ---------------- Grid g(r) computation ------------------

    std::vector<double> compute_pair_correlation_values_once_grid(
        std::size_t bins,
        double maximum_distance,
        double grid_cell_size
    ) const;

    // ---------------- Utility methods ------------------------

    Vector3d apply_minimum_image(Vector3d d, const Domain& domain) const;

    double sample_distance(
        std::size_t number_of_particles,
        const std::vector<Vector3d>& positions,
        std::mt19937_64& rng,
        const Domain& domain
    ) const;

};
