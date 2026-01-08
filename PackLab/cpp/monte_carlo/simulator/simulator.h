#pragma once

#include <cstddef>    // for std::size_t
#include <cstdint>    // for std::uint64_t
#include <vector>     // for std::vector
#include <cmath>      // for std::max
#include <memory>     // for std::shared_ptr, std::unique_ptr

#include "monte_carlo/radius_sampler/radius_sampler.h"
#include "monte_carlo/utils/utils.h"
#include "monte_carlo/statistics/statistics.h"
#include "monte_carlo/simulator/sphere_configuration.h"
#include "monte_carlo/simulator/spatial_grid_index.h"
#include "monte_carlo/result/result.h"
#include "monte_carlo/domain/domain.h"

class Simulator {
public:
    std::shared_ptr<Domain> domain;
    std::shared_ptr<SphereConfiguration> sphere_configuration;
    mutable Statistics statistics;

    Simulator(
        std::shared_ptr<Domain> domain,
        std::shared_ptr<RadiusSampler> radius_sampler,
        std::shared_ptr<Options> options
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
    Attempt to insert a single sphere with a specified radius into the simulation.
    @param radius The radius of the sphere to insert.
    @return True if the insertion was successful, false otherwise.
    */
    bool attempt_single_insertion_with_radius(double radius);


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

    std::shared_ptr<RadiusSampler> radius_sampler;
    std::shared_ptr<Options> options;

    std::mt19937_64 random_generator;

    double maximum_radius_observed = 0.0;
    bool spatial_grid_initialized = false;

    std::unique_ptr<SpatialGridIndex> spatial_grid_index;
};
