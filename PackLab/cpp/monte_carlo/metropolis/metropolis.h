#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>

#include "monte_carlo/domain/domain.h"
#include "monte_carlo/result/result.h"
#include "monte_carlo/simulator/sphere_configuration.h"

struct MetropolisOptions {
    /// Seed for reproducible displacement proposals; zero selects a random seed.
    std::uint64_t random_seed = 0;

    /// A sweep contains one displacement proposal per particle.
    std::size_t number_of_sweeps = 1'000;

    /// Maximum absolute displacement along each Cartesian axis, in metres.
    double maximum_displacement = 0.0;
};

struct MetropolisStatistics {
    std::size_t attempted_moves = 0;
    std::size_t accepted_moves = 0;
    std::size_t rejected_moves = 0;
    std::size_t completed_sweeps = 0;
    double total_runtime_seconds = 0.0;

    [[nodiscard]] double acceptance_rate() const;
};

/*
Equilibrium hard-sphere sampler based on single-particle Metropolis moves.

The particle count, radii, and class labels are fixed.  A trial move is
accepted precisely when it keeps every pair of hard spheres non-overlapping;
this is the Metropolis rule for a hard-sphere system with symmetric proposals.
*/
class MetropolisSimulator {
public:
    std::shared_ptr<MCDomain> domain;
    std::shared_ptr<SphereConfiguration> sphere_configuration;
    MetropolisStatistics statistics;

    MetropolisSimulator(
        std::shared_ptr<MCDomain> domain,
        std::shared_ptr<SphereConfiguration> initial_configuration,
        std::shared_ptr<MetropolisOptions> options
    );

    /// Restore the supplied initial configuration and clear move statistics.
    void reset();

    /// Run the configured number of sweeps and return the equilibrated state.
    Result run();

private:
    std::shared_ptr<SphereConfiguration> initial_configuration;
    std::shared_ptr<MetropolisOptions> options;
    std::mt19937_64 random_generator;

    void validate_initial_configuration() const;
    [[nodiscard]] bool is_inside_nonperiodic_domain(const Vector3d& position, double radius) const;
    [[nodiscard]] bool overlaps_another_sphere(std::size_t moved_index, const Vector3d& candidate) const;
    [[nodiscard]] double squared_center_distance(const Vector3d& a, const Vector3d& b) const;
    [[nodiscard]] Statistics result_statistics() const;
};
