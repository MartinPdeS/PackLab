#include "metropolis.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace {

std::shared_ptr<SphereConfiguration> copy_configuration(const std::shared_ptr<SphereConfiguration>& source) {
    if (!source) {
        throw std::invalid_argument("initial_configuration must be valid.");
    }
    return std::make_shared<SphereConfiguration>(*source);
}

}  // namespace

double MetropolisStatistics::acceptance_rate() const {
    if (attempted_moves == 0) {
        return 0.0;
    }
    return static_cast<double>(accepted_moves) / static_cast<double>(attempted_moves);
}

MetropolisSimulator::MetropolisSimulator(
    std::shared_ptr<MCDomain> input_domain,
    std::shared_ptr<SphereConfiguration> input_configuration,
    std::shared_ptr<MetropolisOptions> input_options
)
    : domain(std::move(input_domain)),
      initial_configuration(copy_configuration(input_configuration)),
      options(std::move(input_options)),
      random_generator(0) {
    if (!domain) {
        throw std::invalid_argument("domain must be valid.");
    }
    if (!options) {
        throw std::invalid_argument("options must be valid.");
    }
    if (!(options->maximum_displacement > 0.0) || !std::isfinite(options->maximum_displacement)) {
        throw std::invalid_argument("maximum_displacement must be finite and positive.");
    }

    validate_initial_configuration();
    reset();
}

void MetropolisSimulator::reset() {
    sphere_configuration = copy_configuration(initial_configuration);
    statistics = MetropolisStatistics{};
    random_generator.seed(options->random_seed == 0 ? std::random_device{}() : options->random_seed);
}

bool MetropolisSimulator::is_inside_nonperiodic_domain(const Vector3d& position, double radius) const {
    if (domain->use_periodic_boundaries) {
        return true;
    }
    return position.x >= radius && position.x <= domain->length_x - radius &&
           position.y >= radius && position.y <= domain->length_y - radius &&
           position.z >= radius && position.z <= domain->length_z - radius;
}

double MetropolisSimulator::squared_center_distance(const Vector3d& a, const Vector3d& b) const {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    if (domain->use_periodic_boundaries) {
        dx = domain->minimum_image_displacement(dx, domain->length_x);
        dy = domain->minimum_image_displacement(dy, domain->length_y);
        dz = domain->minimum_image_displacement(dz, domain->length_z);
    }
    return dx * dx + dy * dy + dz * dz;
}

bool MetropolisSimulator::overlaps_another_sphere(std::size_t moved_index, const Vector3d& candidate) const {
    const double moved_radius = sphere_configuration->radii_values[moved_index];
    for (std::size_t index = 0; index < sphere_configuration->center_positions.size(); ++index) {
        if (index == moved_index) {
            continue;
        }
        const double minimum_distance = moved_radius + sphere_configuration->radii_values[index];
        if (squared_center_distance(candidate, sphere_configuration->center_positions[index]) <
            minimum_distance * minimum_distance) {
            return true;
        }
    }
    return false;
}

void MetropolisSimulator::validate_initial_configuration() const {
    const auto& positions = initial_configuration->center_positions;
    const auto& radii = initial_configuration->radii_values;
    const auto& classes = initial_configuration->class_index_values;

    if (positions.size() != radii.size()) {
        throw std::invalid_argument("initial configuration positions and radii must have the same length.");
    }
    if (!classes.empty() && classes.size() != positions.size()) {
        throw std::invalid_argument("initial configuration class labels must match the number of positions.");
    }

    for (std::size_t index = 0; index < positions.size(); ++index) {
        if (!(radii[index] > 0.0) || !std::isfinite(radii[index])) {
            throw std::invalid_argument("initial configuration radii must be finite and positive.");
        }
        const auto& position = positions[index];
        if (!std::isfinite(position.x) || !std::isfinite(position.y) || !std::isfinite(position.z)) {
            throw std::invalid_argument("initial configuration positions must be finite.");
        }
        if (!is_inside_nonperiodic_domain(position, radii[index])) {
            throw std::invalid_argument("initial configuration contains a sphere outside the non-periodic domain.");
        }
    }

    for (std::size_t first = 0; first < positions.size(); ++first) {
        for (std::size_t second = first + 1; second < positions.size(); ++second) {
            const double minimum_distance = radii[first] + radii[second];
            if (squared_center_distance(positions[first], positions[second]) < minimum_distance * minimum_distance) {
                throw std::invalid_argument("initial configuration contains overlapping spheres.");
            }
        }
    }
}

Statistics MetropolisSimulator::result_statistics() const {
    Statistics output;
    const auto& radii = sphere_configuration->radii_values;
    output.sphere_count = radii.size();
    output.packing_fraction_geometry = sphere_configuration->total_sphere_volume() / domain->volume;
    output.packing_fraction_simulator = output.packing_fraction_geometry;
    if (radii.empty()) {
        return output;
    }

    const auto [minimum, maximum] = std::minmax_element(radii.begin(), radii.end());
    output.radius_min = *minimum;
    output.radius_max = *maximum;
    double sum = 0.0;
    for (const double radius : radii) {
        sum += radius;
    }
    output.radius_mean = sum / static_cast<double>(radii.size());

    std::vector<double> sorted_radii = radii;
    std::sort(sorted_radii.begin(), sorted_radii.end());
    const std::size_t midpoint = sorted_radii.size() / 2;
    output.radius_median = sorted_radii.size() % 2 == 0
        ? 0.5 * (sorted_radii[midpoint - 1] + sorted_radii[midpoint])
        : sorted_radii[midpoint];

    double squared_difference_sum = 0.0;
    for (const double radius : radii) {
        const double difference = radius - output.radius_mean;
        squared_difference_sum += difference * difference;
    }
    output.radius_std = std::sqrt(squared_difference_sum / static_cast<double>(radii.size()));
    return output;
}

Result MetropolisSimulator::run() {
    const auto start = std::chrono::steady_clock::now();
    const std::size_t particle_count = sphere_configuration->center_positions.size();
    std::uniform_real_distribution<double> displacement(-options->maximum_displacement, options->maximum_displacement);

    for (std::size_t sweep = 0; sweep < options->number_of_sweeps; ++sweep) {
        for (std::size_t particle = 0; particle < particle_count; ++particle) {
            ++statistics.attempted_moves;
            const Vector3d previous = sphere_configuration->center_positions[particle];
            Vector3d candidate{
                previous.x + displacement(random_generator),
                previous.y + displacement(random_generator),
                previous.z + displacement(random_generator),
            };

            if (domain->use_periodic_boundaries) {
                candidate = domain->wrap_position_if_periodic(candidate);
            }

            const bool reject = !is_inside_nonperiodic_domain(candidate, sphere_configuration->radii_values[particle]) ||
                                overlaps_another_sphere(particle, candidate);
            if (reject) {
                ++statistics.rejected_moves;
                continue;
            }
            sphere_configuration->center_positions[particle] = candidate;
            ++statistics.accepted_moves;
        }
        ++statistics.completed_sweeps;
    }

    const auto finish = std::chrono::steady_clock::now();
    statistics.total_runtime_seconds += std::chrono::duration<double>(finish - start).count();

    std::size_t number_of_classes = 1;
    for (const int class_index : sphere_configuration->class_index_values) {
        if (class_index < 0) {
            throw std::runtime_error("initial configuration contains a negative class index.");
        }
        number_of_classes = std::max(number_of_classes, static_cast<std::size_t>(class_index + 1));
    }
    return Result(sphere_configuration, domain, result_statistics(), number_of_classes);
}
