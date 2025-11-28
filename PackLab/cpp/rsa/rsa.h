#pragma once

#include <cstddef>    // for std::size_t
#include <cstdint>    // for std::uint64_t
#include <vector>     // for std::vector
#include <cmath>      // for std::max
#include <memory>   // for std::shared_ptr, std::unique_ptr

#include "../radius_sampler/radius_sampler.h"
#include "../utils/utils.h"
#include "../statistics/statistics.h"
#include "sphere_configuration.h"
#include "spatial_grid_index.h"
#include "../result/result.h"
#include "../domain/domain.h"

class Simulator {
public:
    Domain domain;

    Simulator(
        Domain domain,
        std::shared_ptr<RadiusSampler> radius_sampler,
        Options options
    );

    /*
    Reset the simulation to its initial state.
    */
    void reset();

    /*
    Run the simulation.
    @return True if the simulation ran successfully, false otherwise.
    */
    Result run();

    /*
    Attempt to insert a single sphere into the simulation.
    @return True if the insertion was successful, false otherwise.
    */
    bool attempt_single_insertion();

    /*
    Get the current sphere configuration.
    @return A reference to the sphere configuration.
    */
    const SphereConfiguration& sphere_configuration() const { return sphere_configuration_value_; }

    /*
    Get the current simulation statistics.
    @return A reference to the statistics structure.
    */
    const Statistics& statistics() const { return statistics_value_; }

    /*
    Collect the attempted insertion positions.
    @return A reference to the vector of attempted positions.
    */
    const std::vector<Vector3d>& attempted_positions() const { return attempted_positions_values_; }

    /*
    Compute the pair correlation function (g(r)) of the current sphere configuration.
    @param bins Number of bins to use for the histogram.
    @param maximum_pairs Maximum number of random pairs to sample for the calculation.
    @param random_seed Seed for the random number generator used to sample pairs.
    @return A pair of vectors: the first contains the bin centers (r values), and
            the second contains the corresponding g(r) values.
    */
    std::pair<std::vector<double>, std::vector<double>> compute_pair_correlation_function(std::size_t bins, std::size_t maximum_pairs, std::uint64_t random_seed) const;

private:
    /*
    Check if a sphere fits inside the domain considering walls (non-periodic boundaries).
    @param center_position The center position of the sphere.
    @param radius The radius of the sphere.
    @return True if the sphere fits inside the domain, false otherwise.
    */
    bool sphere_fits_inside_domain_if_walls(const Vector3d& center_position, double radius) const;

    /*
    Compute the squared distance between two points considering periodic boundaries.
    @param a The first position.
    @param b The second position.
    @return The squared distance between the two points.
    */
    double periodic_center_distance_squared(const Vector3d& a, const Vector3d& b) const;

    /*
    Compute the squared distance between two points without considering periodic boundaries.
    @param a The first position.
    @param b The second position.
    @return The squared distance between the two points.
    */
    double nonperiodic_center_distance_squared(const Vector3d& a, const Vector3d& b) const;

    /*
    Check if a proposed sphere overlaps with any existing spheres.
    @param center_position The center position of the proposed sphere.
    @param radius The radius of the proposed sphere.
    @return True if there is an overlap, false otherwise.
    */
    bool overlaps_any_existing_sphere(const Vector3d& center_position, double radius) const;

    /*
    Rebuild the spatial grid index if needed based on the new maximum radius.
    @param new_maximum_radius The new maximum radius to consider for grid cell size.
    */
    void rebuild_grid_if_needed(double new_maximum_radius);

private:

    std::shared_ptr<RadiusSampler> radius_sampler_;
    Options options_;

    std::mt19937_64 random_generator_;
    std::chrono::high_resolution_clock::time_point start_time_point_;

    SphereConfiguration sphere_configuration_value_;
    mutable Statistics statistics_value_;

    double maximum_radius_observed_ = 0.0;
    bool spatial_grid_initialized_ = false;

    std::unique_ptr<SpatialGridIndex> spatial_grid_index_;

    std::vector<Vector3d> attempted_positions_values_;
};
