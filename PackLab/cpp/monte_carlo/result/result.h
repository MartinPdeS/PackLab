#pragma once

#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <memory>

#include "monte_carlo/utils/utils.h"
#include "monte_carlo/statistics/statistics.h"
#include "monte_carlo/domain/domain.h"
#include "monte_carlo/simulator/spatial_grid_index.h"
#include "monte_carlo/simulator/sphere_configuration.h"


struct RadialGrid {
    double r_max = 0.0;
    double dr = 0.0;
    std::vector<double> centers;
    std::vector<double> shell_volumes;
};

/*
Result stores a completed RSA configuration along with utilities for
computing total and partial pair correlation functions g(r) and g_ij(r).
*/
class Result {
public:
    std::shared_ptr<Domain> domain;
    std::shared_ptr<SphereConfiguration> sphere_configuration;
    std::size_t number_of_classes;
    Statistics statistics;

    std::vector<double> partial_volume_fractions;
    std::vector<double> partial_volumes;

    void compute_partial_sphere_volumes();


    /*
    Overload including class assignments for polydisperse systems.

    @param sphere_configuration Sphere configuration containing positions, radii, and classes.
    @param domain Simulation domain
    @param number_of_classes Total number of distinct classes.
    */
    Result(const std::shared_ptr<SphereConfiguration>& sphere_configuration, const std::shared_ptr<Domain>& domain, Statistics statistics, std::size_t _number_of_classes)
        :   domain(domain),
            sphere_configuration(sphere_configuration),
            number_of_classes(_number_of_classes),
            statistics(std::move(statistics))
    {
        if (this->sphere_configuration->center_positions.size() != this->sphere_configuration->radii_values.size())
            throw std::invalid_argument("positions and radii must have same length.");

        if (!this->sphere_configuration->class_index_values.empty() && this->sphere_configuration->class_index_values.size() != this->sphere_configuration->center_positions.size())
            throw std::invalid_argument("classes must have same length as positions.");

        this->compute_partial_sphere_volumes();

    }

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

    /*
    Compute partial pair correlation g_ij(r) for class-labeled particles.

    @param n_bins Number of r-bins.
    @param maximum_pairs Maximum number of random pairs to sample.
    @param random_seed RNG seed.

    @return (bin_centers, g_ij[r][i][j]):
            - bin_centers[k] is center of bin k,
            - g_ij[k][i][j] is the partial pair-correlation between class i and j at bin k.
    */
    std::tuple<std::vector<double>, std::vector<std::vector<std::vector<double>>>>
    compute_partial_pair_correlation_function(
        std::size_t n_bins,
        std::size_t maximum_pairs
    ) const;

    /*
        Compute all pair distances along with their class indices.
        @param maximum_pairs Maximum number of random pairs to sample.
        @return (distances, class_i_indices, class_j_indices)
    */
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
    compute_partial_pair_distances(
        std::size_t maximum_pairs
    ) const;

private:

    RadialGrid get_radial_grid(std::size_t number_of_distance_bins) const;

    void validate_partial_pair_inputs(std::size_t number_of_distance_bins) const;

    std::vector<std::size_t> count_particles_per_class() const;

    std::vector<std::vector<std::vector<double>>> make_zero_matrix_double(std::size_t K, std::size_t B) const;

    std::vector<std::vector<std::vector<double>>> make_zero_matrix_expected(std::size_t K, std::size_t B) const;

    double compute_pair_sampling_scale(
        std::size_t N,
        std::size_t examined_unordered_pairs
    ) const;

    std::vector<std::vector<std::vector<std::size_t>>> build_partial_histogram_ordered(
        std::size_t number_of_distance_bins,
        double r_max,
        double dr,
        std::size_t maximum_pairs,
        std::size_t& examined_unordered_pairs
    ) const;

    std::vector<std::vector<std::vector<double>>> compute_g_from_histogram_and_expected(
        const std::vector<std::vector<std::vector<std::size_t>>>& histogram_ordered,
        const std::vector<std::vector<std::vector<double>>>& expected_ordered_counts,
        double sampling_scale
    ) const;

    std::tuple<std::vector<double>, std::vector<std::vector<std::vector<double>>>>
    get_uncorrelated_Cij(
        std::size_t n_bins
    ) const;


    // ---------------- Utility methods ------------------------

    Vector3d apply_minimum_image(Vector3d d) const;

    double sample_distance(
        std::size_t number_of_particles,
        const std::vector<Vector3d>& positions,
        std::mt19937_64& rng
    ) const;


};
