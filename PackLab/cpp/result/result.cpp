#include "result.h"


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

double Result::sample_distance(std::size_t number_of_particles, const std::vector<Vector3d>& positions, std::mt19937_64& random_generator, const Domain& domain)
{
    std::uniform_int_distribution<std::size_t> dist(0, number_of_particles - 1);

    std::size_t i, j;
    do {
        i = dist(random_generator);
        j = dist(random_generator);
    } while (i == j);

    Vector3d d = positions[j] - positions[i];
    if (domain.use_periodic_boundaries()) {
        d = apply_minimum_image(d, domain);
    }

    return d.norm();
}



void Result::compute_pair_correlation_function(std::size_t bins, std::size_t maximum_number_of_pairs, std::uint64_t random_seed) {
    const std::size_t number_of_particles = particle_positions_.size();

    pair_correlation_centers_.clear();
    pair_correlation_values_.clear();

    if (number_of_particles < 2)
        return;

    const double number_density = static_cast<double>(number_of_particles) / domain.volume();

    const double distance_maximum = 0.5 * std::min({domain.length_x(), domain.length_y(), domain.length_z()});
    const double dr = distance_maximum / static_cast<double>(bins);

    std::vector<std::size_t> histogram(bins, 0);

    std::mt19937_64 random_generator(random_seed);

    // sample random pairs
    for (std::size_t k = 0; k < maximum_number_of_pairs; ++k) {
        double r = sample_distance(number_of_particles, particle_positions_, random_generator, domain);

        if (r < distance_maximum) {
            std::size_t index = static_cast<std::size_t>(r / dr);
            if (index < bins)
                histogram[index] += 1;
        }
    }

    // build output arrays
    pair_correlation_centers_.resize(bins);
    pair_correlation_values_.resize(bins);

    const double total_pairs = static_cast<double>(maximum_number_of_pairs);
    const double normalization_factor = (number_of_particles * (number_of_particles - 1) / 2.0) / total_pairs;

    for (std::size_t i = 0; i < bins; ++i) {
        const double radius_inner = i * dr;
        const double radius_outer = (i + 1) * dr;

        const double shell_volume = (4.0 * M_PI / 3.0) * (std::pow(radius_outer, 3) - std::pow(radius_inner, 3));

        pair_correlation_centers_[i] = (radius_inner + radius_outer) * 0.5;
        const double expected = number_density * shell_volume * normalization_factor;

        pair_correlation_values_[i] =
            expected > 0.0 ? histogram[i] / expected : 0.0;
    }
}

