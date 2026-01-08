// simulator.cpp
#include "simulator.h"
#include <iostream>

Simulator::Simulator(std::shared_ptr<Domain> _domain, std::shared_ptr<RadiusSampler> _radius_sampler, std::shared_ptr<Options> _options):
    domain(std::move(_domain)),
    radius_sampler(std::move(_radius_sampler)),
    options(std::move(_options)),
    random_generator(options->random_seed == 0 ? std::random_device{}() : options->random_seed),
    spatial_grid_index(nullptr)
{
    if (!radius_sampler) {
        throw std::invalid_argument("radius_sampler must be valid.");
    }
    reset();
}

void Simulator::reset() {
    // sphere_configuration = SphereConfiguration{};
    sphere_configuration = std::make_shared<SphereConfiguration>();
    this->statistics = Statistics{};

    maximum_radius_observed = 0.0;
    spatial_grid_initialized = false;
    spatial_grid_index.reset();

    this->statistics.start_benchmark();
}

bool Simulator::sphere_fits_inside_domain_if_walls(const Vector3d& center_position, double radius) const {
    if (domain->use_periodic_boundaries)
        return true;

    const double padding = std::max(0.0, options->containment_padding);
    const double effective_radius = radius + padding;
    return (center_position.x >= effective_radius) &&
           (center_position.y >= effective_radius) &&
           (center_position.z >= effective_radius) &&
           (center_position.x <= domain->length_x - effective_radius) &&
           (center_position.y <= domain->length_y - effective_radius) &&
           (center_position.z <= domain->length_z - effective_radius);
}

double Simulator::periodic_center_distance_squared(const Vector3d& a, const Vector3d& b) const {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;

    dx = domain->minimum_image_displacement(dx, domain->length_x);
    dy = domain->minimum_image_displacement(dy, domain->length_y);
    dz = domain->minimum_image_displacement(dz, domain->length_z);

    return dx * dx + dy * dy + dz * dz;
}

double Simulator::nonperiodic_center_distance_squared(const Vector3d& a, const Vector3d& b) const {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

void Simulator::rebuild_grid_if_needed(double new_maximum_radius) {
    if (options->spatial_grid_cell_size > 0.0) {
        if (!spatial_grid_initialized) {

            spatial_grid_index = std::make_unique<SpatialGridIndex>(options->spatial_grid_cell_size, domain);

            spatial_grid_initialized = true;

            for (std::size_t sphere_index = 0; sphere_index < sphere_configuration->center_positions.size(); ++sphere_index)
                spatial_grid_index->insert_sphere(sphere_index, sphere_configuration->center_positions[sphere_index]);

        }
        return;
    }

    const double cell_size_candidate = std::max(1e-12, 2.0 * new_maximum_radius);

    if (!spatial_grid_initialized) {
        spatial_grid_index = std::make_unique<SpatialGridIndex>(cell_size_candidate, domain);
        spatial_grid_initialized = true;
    } else {
        const double previous_size = spatial_grid_index->cell_size();

        if (cell_size_candidate > previous_size * 1.25) {
            spatial_grid_index = std::make_unique<SpatialGridIndex>(cell_size_candidate, domain);
            spatial_grid_index->clear();
            for (std::size_t sphere_index = 0; sphere_index < sphere_configuration->center_positions.size(); ++sphere_index) {
                spatial_grid_index->insert_sphere(sphere_index, sphere_configuration->center_positions[sphere_index]);
            }
        }
    }
}

bool Simulator::overlaps_any_existing_sphere(const Vector3d& center_position, double radius) const {
    if (sphere_configuration->center_positions.empty())
        return false;

    const double extra_separation = std::max(0.0, options->minimum_center_separation_addition);

    const auto neighbor_indices = spatial_grid_index->query_neighbor_sphere_indices(center_position);

    for (std::size_t neighbor_index : neighbor_indices) {
        const Vector3d& neighbor_position = sphere_configuration->center_positions[neighbor_index];
        const double neighbor_radius = sphere_configuration->radii_values[neighbor_index];

        const double required_distance = (radius + neighbor_radius + extra_separation);
        const double required_distance_squared = required_distance * required_distance;

        const double actual_distance_squared = domain->use_periodic_boundaries
            ? periodic_center_distance_squared(center_position, neighbor_position)
            : nonperiodic_center_distance_squared(center_position, neighbor_position);

        if (actual_distance_squared < required_distance_squared) {
            return true;
        }
    }

    return false;
}


bool Simulator::attempt_single_insertion() {
    this->statistics.attempted_insertions += 1;

    const double radius = radius_sampler->sample_radius(random_generator);

    maximum_radius_observed = std::max(maximum_radius_observed, radius);
    rebuild_grid_if_needed(std::max(maximum_radius_observed, radius_sampler->maximum_possible_radius()));

    const double margin = radius + std::max(0.0, options->containment_padding);
    Vector3d proposed_center = domain->sample_uniform_position(random_generator, margin);

    if (!sphere_fits_inside_domain_if_walls(proposed_center, radius)) {
        this->statistics.rejected_insertions += 1;
        this->statistics.consecutive_rejections += 1;
        return false;
    }

    if (overlaps_any_existing_sphere(proposed_center, radius)) {
        this->statistics.rejected_insertions += 1;
        this->statistics.consecutive_rejections += 1;
        return false;
    }

    proposed_center = domain->wrap_position_if_periodic(proposed_center);

    const std::size_t new_index = sphere_configuration->center_positions.size();
    sphere_configuration->center_positions.push_back(proposed_center);

    sphere_configuration->radii_values.push_back(radius);
    const int cls = radius_sampler->bin_index(radius);
    sphere_configuration->class_index_values.push_back(cls);

    spatial_grid_index->insert_sphere(new_index, proposed_center);


    this->statistics.accepted_insertions += 1;
    this->statistics.consecutive_rejections = 0;

    this->statistics.sphere_count = sphere_configuration->radii_values.size();

    this->statistics.packing_fraction_simulator = sphere_configuration->total_sphere_volume() / domain->volume;

    // update radius stats incrementally
    const auto& radii = sphere_configuration->radii_values;
    const double new_radius = radii.back();

    if (this->statistics.sphere_count == 1) {
        this->statistics.radius_min = new_radius;
        this->statistics.radius_max = new_radius;
        this->statistics.radius_mean = new_radius;
    } else {
        this->statistics.radius_min = std::min(this->statistics.radius_min, new_radius);
        this->statistics.radius_max = std::max(this->statistics.radius_max, new_radius);

        // running mean
        const double n = static_cast<double>(this->statistics.sphere_count);
        this->statistics.radius_mean =
            this->statistics.radius_mean * (n - 1.0) / n + new_radius / n;
    }

    return true;
}


bool Simulator::attempt_single_insertion_with_radius(double radius) {
    this->statistics.attempted_insertions += 1;

    maximum_radius_observed = std::max(maximum_radius_observed, radius);
    rebuild_grid_if_needed(std::max(maximum_radius_observed, radius_sampler->maximum_possible_radius()));

    const double margin = radius + std::max(0.0, options->containment_padding);
    Vector3d proposed_center = domain->sample_uniform_position(random_generator, margin);

    if (!sphere_fits_inside_domain_if_walls(proposed_center, radius)) {
        this->statistics.rejected_insertions += 1;
        this->statistics.consecutive_rejections += 1;
        return false;
    }

    if (overlaps_any_existing_sphere(proposed_center, radius)) {
        this->statistics.rejected_insertions += 1;
        this->statistics.consecutive_rejections += 1;
        return false;
    }

    proposed_center = domain->wrap_position_if_periodic(proposed_center);

    const std::size_t new_index = sphere_configuration->center_positions.size();
    sphere_configuration->center_positions.push_back(proposed_center);
    sphere_configuration->radii_values.push_back(radius);

    const int cls = radius_sampler->bin_index(radius);
    sphere_configuration->class_index_values.push_back(cls);

    spatial_grid_index->insert_sphere(new_index, proposed_center);

    this->statistics.accepted_insertions += 1;
    this->statistics.consecutive_rejections = 0;

    this->statistics.sphere_count = sphere_configuration->radii_values.size();
    this->statistics.packing_fraction_simulator =
        sphere_configuration->total_sphere_volume() / domain->volume;

    return true;
}

Result Simulator::run() {
    bool placed_any = false;

    for (std::size_t attempt_index = 0; attempt_index < options->maximum_attempts; ++attempt_index) {
        if (options->maximum_spheres > 0 && sphere_configuration->radii_values.size() >= options->maximum_spheres)
            break;

        if (options->target_packing_fraction > 0.0 && this->statistics.packing_fraction_simulator >= options->target_packing_fraction)
            break;

        if (options->maximum_consecutive_rejections > 0 && this->statistics.consecutive_rejections >= options->maximum_consecutive_rejections)
            break;

        if (!options->enforce_radii_distribution) {
            const bool accepted = this->attempt_single_insertion();
            placed_any = placed_any || accepted;
            continue;
        }

        const double radius = radius_sampler->sample_radius(random_generator);

        bool accepted = false;
        while (!accepted) {
            if (options->maximum_spheres > 0 && sphere_configuration->radii_values.size() >= options->maximum_spheres)
                break;

            if (options->target_packing_fraction > 0.0 && this->statistics.packing_fraction_simulator >= options->target_packing_fraction)
                break;

            if (options->maximum_consecutive_rejections > 0 && this->statistics.consecutive_rejections >= options->maximum_consecutive_rejections)
                break;

            accepted = this->attempt_single_insertion_with_radius(radius);
            placed_any = placed_any || accepted;
        }
    }

    const auto& radii = sphere_configuration->radii_values;
    if (!radii.empty()) {
        std::vector<double> radii_copy = radii;
        std::sort(radii_copy.begin(), radii_copy.end());
        const std::size_t n = radii_copy.size();
        if (n % 2 == 0)
            this->statistics.radius_median = 0.5 * (radii_copy[n / 2 - 1] + radii_copy[n / 2]);
        else
            this->statistics.radius_median = radii_copy[n / 2];

        double mean = this->statistics.radius_mean;
        double variance_sum = 0.0;
        for (double r : radii) {
            const double diff = r - mean;
            variance_sum += diff * diff;
        }
        this->statistics.radius_std = std::sqrt(variance_sum / n);

        this->statistics.packing_fraction_geometry =
            sphere_configuration->total_sphere_volume() / this->domain->volume;
    }

    this->statistics.end_benchmark();

    Result result{
        this->sphere_configuration,
        this->domain,
        this->statistics,
        this->radius_sampler->number_of_bins()
    };

    return result;
}