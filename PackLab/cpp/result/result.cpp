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

double Result::sample_distance(
    std::size_t number_of_particles,
    const std::vector<Vector3d>& positions,
    std::mt19937_64& random_generator,
    const Domain& domain)
{
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
    std::size_t maximum_number_of_pairs,
    std::uint64_t random_seed,
    double maximum_distance)
{
    const std::size_t number_of_particles = particle_positions_.size();

    pair_correlation_centers_.clear();
    pair_correlation_values_.clear();

    if (number_of_particles < 2) return;
    if (bins == 0) return;
    if (maximum_number_of_pairs == 0) return;

    const double Lx = domain.length_x();
    const double Ly = domain.length_y();
    const double Lz = domain.length_z();

    const double safe_distance_maximum = 0.5 * std::min({Lx, Ly, Lz});

    double distance_maximum = maximum_distance;
    if (distance_maximum <= 0.0) {
        distance_maximum = safe_distance_maximum;
    } else {
        distance_maximum = std::min(distance_maximum, safe_distance_maximum);
    }

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

    pair_correlation_centers_.resize(bins);
    pair_correlation_values_.resize(bins);

    const double M = static_cast<double>(maximum_number_of_pairs);

    for (std::size_t i = 0; i < bins; ++i) {
        const double radius_inner = static_cast<double>(i) * dr;
        const double radius_outer = static_cast<double>(i + 1) * dr;

        const double shell_volume =
            (4.0 * PI / 3.0) *
            (std::pow(radius_outer, 3) - std::pow(radius_inner, 3));

        pair_correlation_centers_[i] = 0.5 * (radius_inner + radius_outer);

        const double expected = M * (shell_volume / domain_volume);

        pair_correlation_values_[i] =
            expected > 0.0 ? static_cast<double>(histogram[i]) / expected : 0.0;
    }
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
    if (maximum_number_of_pairs == 0) return values;

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


void Result::compute_pair_correlation_function_mean_and_std(
    std::size_t bins,
    std::size_t maximum_number_of_pairs,
    std::size_t repeats,
    std::uint64_t random_seed,
    double maximum_distance)
{
    const std::size_t number_of_particles = particle_positions_.size();

    pair_correlation_centers_.clear();
    pair_correlation_values_.clear();
    pair_correlation_mean_values_.clear();
    pair_correlation_std_values_.clear();

    if (number_of_particles < 2) return;
    if (bins == 0) return;
    if (maximum_number_of_pairs == 0) return;
    if (repeats == 0) return;

    const double Lx = domain.length_x();
    const double Ly = domain.length_y();
    const double Lz = domain.length_z();

    const double safe_distance_maximum = 0.5 * std::min({Lx, Ly, Lz});

    double distance_maximum = maximum_distance;
    if (distance_maximum <= 0.0) distance_maximum = safe_distance_maximum;
    if (distance_maximum > safe_distance_maximum) distance_maximum = safe_distance_maximum;

    const double dr = distance_maximum / static_cast<double>(bins);

    pair_correlation_centers_.resize(bins);
    for (std::size_t i = 0; i < bins; ++i) {
        const double radius_inner = static_cast<double>(i) * dr;
        const double radius_outer = static_cast<double>(i + 1) * dr;
        pair_correlation_centers_[i] = 0.5 * (radius_inner + radius_outer);
    }

    std::vector<double> mean(bins, 0.0);
    std::vector<double> m2(bins, 0.0);

    auto splitmix64 = [](std::uint64_t& x) -> std::uint64_t {
        x += 0x9e3779b97f4a7c15ULL;
        std::uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    };

    std::uint64_t seed_state = random_seed ? random_seed : 1ULL;

    #pragma omp parallel for
    for (std::size_t r = 0; r < repeats; ++r) {
        const std::uint64_t run_seed = splitmix64(seed_state);

        const std::vector<double> values = compute_pair_correlation_values_once(
            bins,
            maximum_number_of_pairs,
            run_seed,
            distance_maximum
        );

        const double n = static_cast<double>(r + 1);

        for (std::size_t i = 0; i < bins; ++i) {
            const double x = values[i];
            const double delta = x - mean[i];
            mean[i] += delta / n;
            const double delta2 = x - mean[i];
            m2[i] += delta * delta2;
        }
    }

    pair_correlation_mean_values_ = mean;
    pair_correlation_std_values_.assign(bins, 0.0);

    if (repeats > 1) {
        const double denom = static_cast<double>(repeats - 1);
        for (std::size_t i = 0; i < bins; ++i) {
            pair_correlation_std_values_[i] = std::sqrt(m2[i] / denom);
        }
    }

    // optionally keep the last run g(r) in pair_correlation_values_
    pair_correlation_values_ = pair_correlation_mean_values_;
}
