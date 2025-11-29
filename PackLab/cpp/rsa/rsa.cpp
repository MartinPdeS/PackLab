// rsa.cpp
#include "rsa.h"



Simulator::Simulator(Domain domain, std::shared_ptr<RadiusSampler> radius_sampler, Options options):
    domain(std::move(domain)),
    radius_sampler_(std::move(radius_sampler)), options_(options),
    random_generator_(options_.random_seed == 0 ? std::random_device{}() : options_.random_seed),
    spatial_grid_index_(nullptr)
{
    if (!radius_sampler_) {
        throw std::invalid_argument("radius_sampler must be valid.");
    }
    reset();
}

void Simulator::reset() {
    sphere_configuration_value_ = SphereConfiguration{};
    this->statistics = Statistics{};

    maximum_radius_observed_ = 0.0;
    spatial_grid_initialized_ = false;
    spatial_grid_index_.reset();

    attempted_positions_values_.clear();

    this->statistics.start_benchmark();
}

bool Simulator::sphere_fits_inside_domain_if_walls(const Vector3d& center_position, double radius) const {
    if (domain.use_periodic_boundaries())
        return true;

    const double padding = std::max(0.0, options_.containment_padding);
    const double effective_radius = radius + padding;
    return (center_position.x >= effective_radius) &&
           (center_position.y >= effective_radius) &&
           (center_position.z >= effective_radius) &&
           (center_position.x <= domain.length_x() - effective_radius) &&
           (center_position.y <= domain.length_y() - effective_radius) &&
           (center_position.z <= domain.length_z() - effective_radius);
}

double Simulator::periodic_center_distance_squared(const Vector3d& a, const Vector3d& b) const {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;

    dx = domain.minimum_image_displacement(dx, domain.length_x());
    dy = domain.minimum_image_displacement(dy, domain.length_y());
    dz = domain.minimum_image_displacement(dz, domain.length_z());

    return dx * dx + dy * dy + dz * dz;
}

double Simulator::nonperiodic_center_distance_squared(const Vector3d& a, const Vector3d& b) const {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

void Simulator::rebuild_grid_if_needed(double new_maximum_radius) {
    if (options_.spatial_grid_cell_size > 0.0) {
        if (!spatial_grid_initialized_) {

            spatial_grid_index_ = std::make_unique<SpatialGridIndex>(options_.spatial_grid_cell_size, domain);

            spatial_grid_initialized_ = true;

            for (std::size_t sphere_index = 0; sphere_index < sphere_configuration_value_.center_positions_values_.size(); ++sphere_index)
                spatial_grid_index_->insert_sphere(sphere_index, sphere_configuration_value_.center_positions_values_[sphere_index]);

        }
        return;
    }

    const double cell_size_candidate = std::max(1e-12, 2.0 * new_maximum_radius);

    if (!spatial_grid_initialized_) {
        spatial_grid_index_ = std::make_unique<SpatialGridIndex>(cell_size_candidate, domain);
        spatial_grid_initialized_ = true;
    } else {
        const double previous_size = spatial_grid_index_->cell_size();

        if (cell_size_candidate > previous_size * 1.25) {
            spatial_grid_index_ = std::make_unique<SpatialGridIndex>(cell_size_candidate, domain);
            spatial_grid_index_->clear();
            for (std::size_t sphere_index = 0; sphere_index < sphere_configuration_value_.center_positions_values_.size(); ++sphere_index) {
                spatial_grid_index_->insert_sphere(sphere_index, sphere_configuration_value_.center_positions_values_[sphere_index]);
            }
        }
    }
}

bool Simulator::overlaps_any_existing_sphere(const Vector3d& center_position, double radius) const {
    if (sphere_configuration_value_.center_positions_values_.empty()) {
        return false;
    }

    const double extra_separation = std::max(0.0, options_.minimum_center_separation_addition);

    const auto neighbor_indices = spatial_grid_index_->query_neighbor_sphere_indices(center_position);

    for (std::size_t neighbor_index : neighbor_indices) {
        const Vector3d& neighbor_position = sphere_configuration_value_.center_positions_values_[neighbor_index];
        const double neighbor_radius = sphere_configuration_value_.radii_values_[neighbor_index];

        const double required_distance = (radius + neighbor_radius + extra_separation);
        const double required_distance_squared = required_distance * required_distance;

        const double actual_distance_squared = domain.use_periodic_boundaries()
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

    const double radius = radius_sampler_->sample_radius(random_generator_);
    maximum_radius_observed_ = std::max(maximum_radius_observed_, radius);
    rebuild_grid_if_needed(std::max(maximum_radius_observed_, radius_sampler_->maximum_possible_radius()));

    const double margin = radius + std::max(0.0, options_.containment_padding);
    Vector3d proposed_center = domain.sample_uniform_position(random_generator_, margin);

    if (options_.store_attempt_positions) {
        attempted_positions_values_.push_back(proposed_center);
    }

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

    proposed_center = domain.wrap_position_if_periodic(proposed_center);

    const std::size_t new_index = sphere_configuration_value_.center_positions_values_.size();
    sphere_configuration_value_.center_positions_values_.push_back(proposed_center);
    sphere_configuration_value_.radii_values_.push_back(radius);

    spatial_grid_index_->insert_sphere(new_index, proposed_center);


    this->statistics.accepted_insertions += 1;
    this->statistics.consecutive_rejections = 0;

    this->statistics.sphere_count = sphere_configuration_value_.radii_values_.size();

    this->statistics.packing_fraction_simulator =
        sphere_configuration_value_.total_sphere_volume() / domain.volume();

    // update radius stats incrementally
    const auto& radii = sphere_configuration_value_.radii_values_;
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

Result Simulator::run() {
    bool placed_any = false;

    for (std::size_t attempt_index = 0; attempt_index < options_.maximum_attempts; ++attempt_index)
    {
        if (options_.maximum_spheres > 0 && sphere_configuration_value_.radii_values_.size() >= options_.maximum_spheres)
            break;

        if (options_.target_packing_fraction > 0.0 && this->statistics.packing_fraction_simulator >= options_.target_packing_fraction)
            break;

        if (options_.maximum_consecutive_rejections > 0 && this->statistics.consecutive_rejections >= options_.maximum_consecutive_rejections)
            break;


        const bool accepted = attempt_single_insertion();
        placed_any = placed_any || accepted;
    }

    // compute final radius statistics that require full array
    const auto& radii = sphere_configuration_value_.radii_values_;
    if (!radii.empty()) {
        // compute median
        std::vector<double> radii_copy = radii;
        std::sort(radii_copy.begin(), radii_copy.end());
        const std::size_t n = radii_copy.size();
        if (n % 2 == 0) {
            this->statistics.radius_median =
                0.5 * (radii_copy[n / 2 - 1] + radii_copy[n / 2]);
        } else {
            this->statistics.radius_median = radii_copy[n / 2];
        }

        // compute standard deviation
        double mean = this->statistics.radius_mean;
        double variance_sum = 0.0;
        for (double r : radii) {
            const double diff = r - mean;
            variance_sum += diff * diff;
        }
        this->statistics.radius_std = std::sqrt(variance_sum / n);

        // set geometric packing fraction only once
        this->statistics.packing_fraction_geometry =
            sphere_configuration_value_.total_sphere_volume() /
            domain.volume();
    }

    // runtime
    this->statistics.end_benchmark();

    Result result{
        sphere_configuration_value_.center_positions_values_,
        sphere_configuration_value_.radii_values_,
        domain
    };

    return result;
}



std::pair<std::vector<double>, std::vector<double>>
Simulator::compute_pair_correlation_function(std::size_t bins, std::size_t maximum_pairs, std::uint64_t random_seed) const
{
    const auto number_of_particles = sphere_configuration_value_.center_positions_values_.size();
    if (number_of_particles < 2) {
        return {{}, {}};
    }

    std::mt19937_64 rng(random_seed);
    std::uniform_int_distribution<std::size_t> distribution(0, number_of_particles - 1);

    const double Lx = domain.length_x();
    const double Ly = domain.length_y();
    const double Lz = domain.length_z();

    const double box_volume = Lx * Ly * Lz;
    const double number_density = number_of_particles / box_volume;

    const double distance_maximum = 0.5 * std::min({Lx, Ly, Lz});

    std::vector<double> distances;
    distances.reserve(maximum_pairs);

    std::size_t pair_count = 0;

    while (pair_count < maximum_pairs) {
        const std::size_t i = distribution(rng);
        const std::size_t j = distribution(rng);
        if (i == j) continue;

        Vector3d d = sphere_configuration_value_.center_positions_values_[j] - sphere_configuration_value_.center_positions_values_[i];

        if (domain.use_periodic_boundaries()) {
            const double hx = 0.5 * Lx;
            const double hy = 0.5 * Ly;
            const double hz = 0.5 * Lz;

            if      (d.x >  hx) d.x -= Lx;
            else if (d.x < -hx) d.x += Lx;

            if      (d.y >  hy) d.y -= Ly;
            else if (d.y < -hy) d.y += Ly;

            if      (d.z >  hz) d.z -= Lz;
            else if (d.z < -hz) d.z += Lz;
        }

        const double r = d.norm();
        if (r < distance_maximum) {
            distances.push_back(r);
            pair_count += 1;
        }
    }

    if (distances.empty()) {
        return {{}, {}};
    }

    std::vector<std::size_t> histogram(bins, 0);
    const double dr = distance_maximum / bins;

    for (double r : distances) {
        std::size_t index = static_cast<std::size_t>(r / dr);
        if (index < bins) histogram[index] += 1;
    }

    std::vector<double> centers(bins);
    std::vector<double> g_r(bins);

    const double number_of_possible_pairs =
        number_of_particles * (number_of_particles - 1.0) * 0.5;

    const double normalization_factor =
        static_cast<double>(pair_count) / number_of_possible_pairs;

    for (std::size_t i = 0; i < bins; ++i) {
        const double r_inner = i * dr;
        const double r_outer = (i + 1) * dr;

        const double shell_volume =
            (4.0 * PI / 3.0) * (std::pow(r_outer, 3) - std::pow(r_inner, 3));

        const double expected =
            number_density * shell_volume * normalization_factor;

        centers[i] = 0.5 * (r_inner + r_outer);
        g_r[i] = expected > 0.0 ? histogram[i] / expected : 0.0;
    }

    return {centers, g_r};
}

