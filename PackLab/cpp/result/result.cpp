#include "result.h"




namespace {
struct PairCorrelationParameters {
    double maximum_distance = 0.0;
    double radial_bin_width = 0.0;
    bool periodic = true;
};

double resolve_grid_cell_size(double grid_cell_size, double maximum_distance, const Domain& domain)
{
    const double length_x = domain.length_x();
    const double length_y = domain.length_y();
    const double length_z = domain.length_z();

    const bool periodic = domain.use_periodic_boundaries();

    const double periodic_safe_maximum = 0.5 * std::min({length_x, length_y, length_z});
    const double nonperiodic_possible_maximum = std::sqrt(length_x * length_x + length_y * length_y + length_z * length_z);

    const double maximum_allowed = periodic ? periodic_safe_maximum : nonperiodic_possible_maximum;

    double resolved = grid_cell_size;

    if (resolved <= 0.0) {
        resolved = maximum_distance;
    } else {
        const double minimum_reasonable = maximum_distance / 4.0;
        if (resolved < minimum_reasonable) resolved = minimum_reasonable;
        if (resolved > maximum_allowed) resolved = maximum_allowed;
    }

    return resolved;
}


} // namespace



void Result::fill_pair_correlation_histogram_monte_carlo(
    std::vector<std::size_t>& histogram,
    std::size_t bins,
    std::size_t maximum_number_of_pairs,
    std::uint64_t random_seed,
    double maximum_distance,
    double radial_bin_width) const
{
    const std::size_t number_of_particles = particle_positions_.size();
    histogram.assign(bins, 0);

    std::mt19937_64 random_generator(random_seed);

    for (std::size_t k = 0; k < maximum_number_of_pairs; ++k) {
        const double r = this->sample_distance(
            number_of_particles,
            particle_positions_,
            random_generator,
            domain
        );

        if (r >= maximum_distance) continue;

        const std::size_t index = static_cast<std::size_t>(r / radial_bin_width);
        if (index < bins) {
            histogram[index] += 1;
        }
    }
}

void Result::fill_pair_correlation_histogram_grid(
    std::vector<std::size_t>& histogram,
    std::size_t bins,
    double maximum_distance,
    double radial_bin_width,
    double grid_cell_size) const
{
    const std::size_t number_of_particles = particle_positions_.size();
    histogram.assign(bins, 0);

    const double resolved_grid_cell_size = resolve_grid_cell_size(grid_cell_size, maximum_distance, domain);

    SpatialGridIndex spatial_grid_index(resolved_grid_cell_size, domain);

    for (std::size_t i = 0; i < number_of_particles; ++i) {
        spatial_grid_index.insert_sphere(i, particle_positions_[i]);
    }

    const std::int64_t cell_range =
        static_cast<std::int64_t>(std::ceil(maximum_distance / spatial_grid_index.cell_size()));

    std::vector<std::size_t> neighbor_indices;
    neighbor_indices.reserve(512);

    const bool periodic = domain.use_periodic_boundaries();

    const double length_x = domain.length_x();
    const double length_y = domain.length_y();
    const double length_z = domain.length_z();

    for (std::size_t i = 0; i < number_of_particles; ++i) {
        neighbor_indices.clear();

        spatial_grid_index.append_neighbor_sphere_indices(
            particle_positions_[i],
            cell_range,
            neighbor_indices
        );

        const Vector3d& pi = particle_positions_[i];

        for (std::size_t j : neighbor_indices) {
            if (j <= i) continue;

            const Vector3d& pj = particle_positions_[j];

            double dx = pj.x - pi.x;
            double dy = pj.y - pi.y;
            double dz = pj.z - pi.z;

            if (periodic) {
                dx = domain.minimum_image_displacement(dx, length_x);
                dy = domain.minimum_image_displacement(dy, length_y);
                dz = domain.minimum_image_displacement(dz, length_z);
            }

            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r >= maximum_distance) continue;

            const std::size_t index = static_cast<std::size_t>(r / radial_bin_width);
            if (index < bins) {
                histogram[index] += 1;
            }
        }
    }
}















Vector3d Result::apply_minimum_image(Vector3d d, const Domain& domain) const {
    const double hx = 0.5 * domain.length_x();
    const double hy = 0.5 * domain.length_y();
    const double hz = 0.5 * domain.length_z();

    if      (d.x >  hx) d.x -= domain.length_x();
    else if (d.x < -hx) d.x += domain.length_x();

    if      (d.y >  hy) d.y -= domain.length_y();
    else if (d.y < -hy) d.y += domain.length_y();

    if      (d.z >  hz) d.z -= domain.length_z();
    else if (d.z < -hz) d.z += domain.length_z();

    return d;
}

double Result::sample_distance(std::size_t number_of_particles, const std::vector<Vector3d>& positions, std::mt19937_64& random_generator, const Domain& domain) const {
    std::uniform_int_distribution<std::size_t> dist(0, number_of_particles - 1);

    std::size_t i = 0;
    std::size_t j = 0;

    do {
        i = dist(random_generator);
        j = dist(random_generator);
    } while (i == j);

    Vector3d displacement = positions[j] - positions[i];

    if (domain.use_periodic_boundaries()) {
        displacement = apply_minimum_image(displacement, domain);
    }

    return displacement.norm();
}

void Result::compute_pair_correlation_function(
    std::size_t bins,
    double maximum_distance,
    const std::string& method,
    std::size_t maximum_number_of_pairs,
    std::uint64_t random_seed,
    double grid_cell_size
)
{
    pair_correlation_centers_.clear();
    pair_correlation_values_.clear();

    const std::size_t number_of_particles = particle_positions_.size();

    if (number_of_particles < 2) return;
    if (bins == 0) return;

    const double length_x = domain.length_x();
    const double length_y = domain.length_y();
    const double length_z = domain.length_z();

    const bool periodic = domain.use_periodic_boundaries();

    const double periodic_safe_maximum = 0.5 * std::min({length_x, length_y, length_z});
    const double nonperiodic_possible_maximum = std::sqrt(length_x * length_x + length_y * length_y + length_z * length_z);

    const double default_maximum_distance = periodic ? periodic_safe_maximum : nonperiodic_possible_maximum;

    double distance_maximum = maximum_distance;

    if (distance_maximum <= 0.0) {
        distance_maximum = default_maximum_distance;
    } else {
        if (periodic && distance_maximum > periodic_safe_maximum) distance_maximum = periodic_safe_maximum;
        if (!periodic && distance_maximum > nonperiodic_possible_maximum) distance_maximum = nonperiodic_possible_maximum;
    }

    const double radial_bin_width = distance_maximum / static_cast<double>(bins);

    pair_correlation_centers_.resize(bins);
    for (std::size_t i = 0; i < bins; ++i) {
        const double radius_inner = static_cast<double>(i) * radial_bin_width;
        const double radius_outer = static_cast<double>(i + 1) * radial_bin_width;
        pair_correlation_centers_[i] = 0.5 * (radius_inner + radius_outer);
    }

    std::string method_normalized = method;
    std::transform(
        method_normalized.begin(),
        method_normalized.end(),
        method_normalized.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); }
    );

    if (method_normalized == "grid") {
        pair_correlation_values_ = compute_pair_correlation_values_once_grid(
            bins,
            distance_maximum,
            grid_cell_size
        );
        return;
    }

    if (method_normalized == "monte-carlo" || method_normalized == "mc") {
        if (maximum_number_of_pairs == 0) {
            throw std::invalid_argument("maximum_number_of_pairs must be positive when method is monte-carlo.");
        }

        pair_correlation_values_ = compute_pair_correlation_values_once(
            bins,
            maximum_number_of_pairs,
            random_seed,
            distance_maximum
        );
        return;
    }

    throw std::invalid_argument("method must be 'grid' or 'monte-carlo'.");
}


std::vector<double> Result::compute_pair_correlation_values_once(
    std::size_t bins,
    std::size_t maximum_number_of_pairs,
    std::uint64_t random_seed,
    double distance_maximum)
{
    const std::size_t number_of_particles = particle_positions_.size();
    std::vector<double> values(bins, 0.0);

    if (number_of_particles < 2) return values;
    if (bins == 0) return values;
    if (maximum_number_of_pairs == 0)
        return values;

    const double domain_volume = domain.volume();
    const double dr = distance_maximum / static_cast<double>(bins);

    std::vector<std::size_t> histogram(bins, 0);

    std::mt19937_64 random_generator(random_seed);

    for (std::size_t k = 0; k < maximum_number_of_pairs; ++k) {
        const double r = sample_distance(
            number_of_particles,
            particle_positions_,
            random_generator,
            domain
        );

        if (r >= distance_maximum) continue;

        const std::size_t index = static_cast<std::size_t>(r / dr);
        if (index < bins) {
            histogram[index] += 1;
        }
    }

    const double M = static_cast<double>(maximum_number_of_pairs);

    for (std::size_t i = 0; i < bins; ++i) {
        const double radius_inner = static_cast<double>(i) * dr;
        const double radius_outer = static_cast<double>(i + 1) * dr;

        const double shell_volume =
            (4.0 * PI / 3.0) *
            (std::pow(radius_outer, 3) - std::pow(radius_inner, 3));

        const double expected = M * (shell_volume / domain_volume);
        values[i] = expected > 0.0 ? static_cast<double>(histogram[i]) / expected : 0.0;
    }

    return values;
}


std::vector<double> Result::compute_pair_correlation_values_once(
    std::size_t bins,
    std::size_t maximum_number_of_pairs,
    std::uint64_t random_seed,
    double maximum_distance
) const
{
    const std::size_t number_of_particles = particle_positions_.size();
    std::vector<double> values(bins, 0.0);

    if (number_of_particles < 2) return values;
    if (bins == 0) return values;
    if (maximum_number_of_pairs == 0) return values;

    const double length_x = domain.length_x();
    const double length_y = domain.length_y();
    const double length_z = domain.length_z();

    const bool periodic = domain.use_periodic_boundaries();

    const double periodic_safe_maximum = 0.5 * std::min({length_x, length_y, length_z});
    const double nonperiodic_possible_maximum = std::sqrt(length_x * length_x + length_y * length_y + length_z * length_z);

    const double default_maximum_distance = periodic ? periodic_safe_maximum : nonperiodic_possible_maximum;

    double distance_maximum = maximum_distance;
    if (distance_maximum <= 0.0) {
        distance_maximum = default_maximum_distance;
    } else {
        if (periodic && distance_maximum > periodic_safe_maximum) distance_maximum = periodic_safe_maximum;
        if (!periodic && distance_maximum > nonperiodic_possible_maximum) distance_maximum = nonperiodic_possible_maximum;
    }

    const double domain_volume = domain.volume();
    const double radial_bin_width = distance_maximum / static_cast<double>(bins);

    std::vector<std::size_t> histogram(bins, 0);

    std::mt19937_64 random_generator(random_seed);

    std::uniform_int_distribution<std::size_t> dist(0, number_of_particles - 1);

    const bool use_minimum_image = periodic;

    for (std::size_t k = 0; k < maximum_number_of_pairs; ++k) {
        std::size_t i = dist(random_generator);
        std::size_t j = dist(random_generator);
        while (j == i) {
            j = dist(random_generator);
        }

        const Vector3d& pi = particle_positions_[i];
        const Vector3d& pj = particle_positions_[j];

        double dx = pj.x - pi.x;
        double dy = pj.y - pi.y;
        double dz = pj.z - pi.z;

        if (use_minimum_image) {
            dx = domain.minimum_image_displacement(dx, length_x);
            dy = domain.minimum_image_displacement(dy, length_y);
            dz = domain.minimum_image_displacement(dz, length_z);
        }

        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);

        if (r >= distance_maximum) continue;

        const std::size_t index = static_cast<std::size_t>(r / radial_bin_width);
        if (index < bins) {
            histogram[index] += 1;
        }
    }

    const double total_pairs_used = static_cast<double>(maximum_number_of_pairs);

    for (std::size_t i = 0; i < bins; ++i) {
        const double radius_inner = static_cast<double>(i) * radial_bin_width;
        const double radius_outer = static_cast<double>(i + 1) * radial_bin_width;

        const double shell_volume =
            (4.0 * PI / 3.0) * (std::pow(radius_outer, 3) - std::pow(radius_inner, 3));

        const double expected = total_pairs_used * (shell_volume / domain_volume);

        values[i] = expected > 0.0 ? static_cast<double>(histogram[i]) / expected : 0.0;
    }

    return values;
}


std::vector<double> Result::compute_pair_correlation_values_once_grid(
    std::size_t bins,
    double maximum_distance,
    double grid_cell_size
) const
{
    const std::size_t number_of_particles = particle_positions_.size();
    std::vector<double> values(bins, 0.0);

    if (number_of_particles < 2) return values;
    if (bins == 0) return values;

    const double length_x = domain.length_x();
    const double length_y = domain.length_y();
    const double length_z = domain.length_z();

    const bool periodic = domain.use_periodic_boundaries();

    const double periodic_safe_maximum = 0.5 * std::min({length_x, length_y, length_z});
    const double nonperiodic_possible_maximum = std::sqrt(length_x * length_x + length_y * length_y + length_z * length_z);

    const double default_maximum_distance = periodic ? periodic_safe_maximum : nonperiodic_possible_maximum;

    double distance_maximum = maximum_distance;

    if (distance_maximum <= 0.0) {
        distance_maximum = default_maximum_distance;
    } else {
        if (periodic && distance_maximum > periodic_safe_maximum) distance_maximum = periodic_safe_maximum;
        if (!periodic && distance_maximum > nonperiodic_possible_maximum) distance_maximum = nonperiodic_possible_maximum;
    }

    double resolved_grid_cell_size = grid_cell_size;

    if (resolved_grid_cell_size <= 0.0) {
        resolved_grid_cell_size = distance_maximum;
    } else {
        const double minimum_reasonable = distance_maximum / 4.0;
        if (resolved_grid_cell_size < minimum_reasonable) resolved_grid_cell_size = minimum_reasonable;

        const double maximum_allowed = periodic ? periodic_safe_maximum : nonperiodic_possible_maximum;
        if (resolved_grid_cell_size > maximum_allowed) resolved_grid_cell_size = maximum_allowed;
    }

    const double domain_volume = domain.volume();
    const double radial_bin_width = distance_maximum / static_cast<double>(bins);

    std::vector<std::size_t> histogram(bins, 0);

    SpatialGridIndex spatial_grid_index(resolved_grid_cell_size, domain);

    for (std::size_t i = 0; i < number_of_particles; ++i) {
        spatial_grid_index.insert_sphere(i, particle_positions_[i]);
    }

    const std::int64_t cell_range =
        static_cast<std::int64_t>(std::ceil(distance_maximum / spatial_grid_index.cell_size()));

    std::vector<std::size_t> neighbor_indices;
    neighbor_indices.reserve(512);

    for (std::size_t i = 0; i < number_of_particles; ++i) {
        neighbor_indices.clear();

        spatial_grid_index.append_neighbor_sphere_indices(
            particle_positions_[i],
            cell_range,
            neighbor_indices
        );

        const Vector3d& pi = particle_positions_[i];

        for (std::size_t j : neighbor_indices) {
            if (j <= i) continue;

            const Vector3d& pj = particle_positions_[j];

            double dx = pj.x - pi.x;
            double dy = pj.y - pi.y;
            double dz = pj.z - pi.z;

            if (periodic) {
                dx = domain.minimum_image_displacement(dx, length_x);
                dy = domain.minimum_image_displacement(dy, length_y);
                dz = domain.minimum_image_displacement(dz, length_z);
            }

            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (r >= distance_maximum) continue;

            const std::size_t index = static_cast<std::size_t>(r / radial_bin_width);
            if (index < bins) {
                histogram[index] += 1;
            }
        }
    }

    const double total_possible_pairs =
        static_cast<double>(number_of_particles) * static_cast<double>(number_of_particles - 1) * 0.5;

    for (std::size_t i = 0; i < bins; ++i) {
        const double radius_inner = static_cast<double>(i) * radial_bin_width;
        const double radius_outer = static_cast<double>(i + 1) * radial_bin_width;

        const double shell_volume =
            (4.0 * PI / 3.0) * (std::pow(radius_outer, 3) - std::pow(radius_inner, 3));

        const double expected = total_possible_pairs * (shell_volume / domain_volume);

        values[i] = expected > 0.0 ? static_cast<double>(histogram[i]) / expected : 0.0;
    }

    return values;
}
