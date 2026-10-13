#include "result.h"
#include "pair_sampling.h"

// -----------------------------------------------------------------------------
// Minimum image helper
// -----------------------------------------------------------------------------
Vector3d Result::apply_minimum_image(Vector3d d) const {
    d.x = domain->minimum_image_displacement(d.x, domain->length_x);
    d.y = domain->minimum_image_displacement(d.y, domain->length_y);
    d.z = domain->minimum_image_displacement(d.z, domain->length_z);

    return d;
}

// -----------------------------------------------------------------------------
// Random pair distance sampler
// -----------------------------------------------------------------------------
double Result::sample_distance(
    std::size_t number_of_particles,
    const std::vector<Vector3d>& positions,
    std::mt19937_64& rng) const
{
    std::uniform_int_distribution<std::size_t> pick(0, number_of_particles - 1);

    std::size_t i = 0, j = 0;
    while ((i = pick(rng)) == (j = pick(rng))) {}

    Vector3d d = positions[j] - positions[i];

    if (domain->use_periodic_boundaries)
        d = apply_minimum_image(d);

    return d.norm();
}

// -----------------------------------------------------------------------------
// Monte Carlo g(r) computation
// -----------------------------------------------------------------------------
void Result::validate_partial_pair_inputs(std::size_t number_of_distance_bins) const
{
    const std::size_t N = this->sphere_configuration->center_positions.size();

    if (N < 2)
        throw std::runtime_error("Need at least two particles to compute partial pair correlation.");

    if (number_of_distance_bins == 0)
        throw std::runtime_error("number_of_distance_bins must be > 0.");

    if (this->sphere_configuration->class_index_values.empty())
        throw std::runtime_error("class_index_values must be set for partial pair correlation.");

    if (this->sphere_configuration->class_index_values.size() != N)
        throw std::runtime_error("class_index_values {" + std::to_string(this->sphere_configuration->class_index_values.size()) + "} and center_positions {" + std::to_string(N) + "} must have the same length.");

    if (number_of_classes <= 0)
        throw std::runtime_error("number_of_classes must be > 0.");

    if (domain->volume <= 0.0)
        throw std::runtime_error("domain->volume must be > 0.");
}

RadialGrid Result::get_radial_grid(std::size_t number_of_distance_bins) const
{
    RadialGrid grid;

    grid.r_max = 0.5 * std::min({domain->length_x, domain->length_y, domain->length_z});
    grid.dr = grid.r_max / static_cast<double>(number_of_distance_bins);

    grid.centers.resize(number_of_distance_bins);
    grid.shell_volumes.resize(number_of_distance_bins);

    for (std::size_t b = 0; b < number_of_distance_bins; ++b) {
        const double r0 = static_cast<double>(b) * grid.dr;
        const double r1 = static_cast<double>(b + 1) * grid.dr;

        grid.centers[b] = (static_cast<double>(b) + 0.5) * grid.dr;

        grid.shell_volumes[b] =
            (4.0 * PI / 3.0) * (std::pow(r1, 3) - std::pow(r0, 3));
    }

    return grid;
}


// -----------------------------------------------------------------------------
// Partial g_ij(r) for class-labeled polydisperse system
// -----------------------------------------------------------------------------
std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
Result::compute_partial_pair_distances(std::size_t maximum_pairs, std::uint64_t random_seed) const
{
    const std::size_t N = this->sphere_configuration->center_positions.size();
    if (N < 2) {
        return { {}, {}, {} };
    }

    if (this->sphere_configuration->class_index_values.empty())
        throw std::runtime_error(
            "compute_partial_pair_distances requires sphere_configuration->class_index_values to be set."
        );

    if (this->sphere_configuration->class_index_values.size() != N)
        throw std::runtime_error(
            "class_index_values and center_positions must have the same length."
        );

    const double r_max = 0.5 * std::min({domain->length_x, domain->length_y, domain->length_z});
    const std::size_t total = N * (N - 1) / 2;
    const std::size_t samples = maximum_pairs == 0 ? total : std::min(maximum_pairs, total);
    std::vector<double> distances;
    std::vector<int> class_i, class_j;
    distances.reserve(samples);
    class_i.reserve(samples);
    class_j.reserve(samples);
    const auto& positions = sphere_configuration->center_positions;
    const auto& classes = sphere_configuration->class_index_values;
    visit_unordered_pairs(N, maximum_pairs, random_seed, [&](std::size_t first, std::size_t second) {
        const auto delta = apply_minimum_image(positions[second] - positions[first]);
        const double distance = delta.norm();
        if (distance >= r_max) return;
        distances.push_back(distance);
        class_i.push_back(classes[first]);
        class_j.push_back(classes[second]);
    });

    return { distances, class_i, class_j };
}


void Result::compute_partial_sphere_volumes()
{
    this->partial_volumes = std::vector<double>(this->number_of_classes, 0.0);

    double prefactor = (4.0 / 3.0) * PI;

    if (this->sphere_configuration->class_index_values.empty()) {
        for (double r : sphere_configuration->radii_values)
            this->partial_volumes[0] += prefactor * r * r * r;
    }



    if (this->sphere_configuration->class_index_values.size() != sphere_configuration->radii_values.size())
        throw std::runtime_error("class_index_values and radii_values must have the same length.");

    for (std::size_t n = 0; n < sphere_configuration->radii_values.size(); ++n) {
        const int c = this->sphere_configuration->class_index_values[n];
        if (c < 0 || static_cast<std::size_t>(c) >= this->number_of_classes) {
            throw std::runtime_error("Invalid class index in class_index_values.");
        }
        this->partial_volumes[c] += prefactor * this->sphere_configuration->radii_values[n] * this->sphere_configuration->radii_values[n] * this->sphere_configuration->radii_values[n];
    }

    for (std::size_t c = 0; c < this->number_of_classes; ++c) {
        this->partial_volume_fractions.push_back(this->partial_volumes[c] / this->domain->volume);
    }
}

std::vector<std::size_t> Result::count_particles_per_class() const
{
    const std::size_t K = static_cast<std::size_t>(number_of_classes);
    std::vector<std::size_t> counts(K, 0);

    for (std::size_t c : this->sphere_configuration->class_index_values) {
        if (c < 0 || c >= number_of_classes)
            throw std::runtime_error("Invalid class index in class_index_values.");
        counts[static_cast<std::size_t>(c)] += 1;
    }

    return counts;
}

std::vector<std::vector<std::vector<double>>>
Result::make_zero_matrix_double(std::size_t K, std::size_t B) const
{
    return std::vector<std::vector<std::vector<double>>>(
        K,
        std::vector<std::vector<double>>(K, std::vector<double>(B, 0.0))
    );
}

std::vector<std::vector<std::vector<double>>>
Result::make_zero_matrix_expected(std::size_t K, std::size_t B) const
{
    return make_zero_matrix_double(K, B);
}

double Result::compute_pair_sampling_scale(std::size_t N, std::size_t examined_unordered_pairs) const
{
    const std::size_t total_unordered_pairs = (N * (N - 1)) / 2;

    if (examined_unordered_pairs == 0)
        return 1.0;

    if (examined_unordered_pairs >= total_unordered_pairs)
        return 1.0;

    return static_cast<double>(total_unordered_pairs) / static_cast<double>(examined_unordered_pairs);
}

std::vector<std::vector<std::vector<std::size_t>>> Result::build_partial_histogram_ordered(
    std::size_t number_of_distance_bins,
    double r_max,
    double dr,
    std::size_t maximum_pairs,
    std::uint64_t random_seed,
    std::size_t& examined_unordered_pairs) const
{
    const std::size_t N = this->sphere_configuration->center_positions.size();
    const std::size_t K = static_cast<std::size_t>(number_of_classes);

    std::vector<std::vector<std::vector<std::size_t>>> histogram_ordered(
        K,
        std::vector<std::vector<std::size_t>>(K, std::vector<std::size_t>(number_of_distance_bins, 0))
    );

    examined_unordered_pairs = 0;
    visit_unordered_pairs(N, maximum_pairs, random_seed, [&](std::size_t first, std::size_t second) {
        ++examined_unordered_pairs;
        const auto ic = static_cast<std::size_t>(sphere_configuration->class_index_values[first]);
        const auto jc = static_cast<std::size_t>(sphere_configuration->class_index_values[second]);
        const auto delta = apply_minimum_image(sphere_configuration->center_positions[second] -
                                               sphere_configuration->center_positions[first]);
        const double distance = delta.norm();
        if (distance >= r_max) return;
        const auto bin = static_cast<std::size_t>(distance / dr);
        if (bin >= number_of_distance_bins) return;
        ++histogram_ordered[ic][jc][bin];
        ++histogram_ordered[jc][ic][bin];
    });

    return histogram_ordered;
}

std::vector<std::vector<std::vector<double>>>
Result::compute_g_from_histogram_and_expected(
    const std::vector<std::vector<std::vector<std::size_t>>>& histogram_ordered,
    const std::vector<std::vector<std::vector<double>>>& expected_ordered_counts,
    double sampling_scale) const
{
    const std::size_t K = histogram_ordered.size();
    const std::size_t B = histogram_ordered[0][0].size();

    std::vector<std::vector<std::vector<double>>> g = make_zero_matrix_double(K, B);

    for (std::size_t i = 0; i < K; ++i) {
        for (std::size_t j = 0; j < K; ++j) {
            for (std::size_t b = 0; b < B; ++b) {
                const double observed =
                    static_cast<double>(histogram_ordered[i][j][b]) * sampling_scale;

                const double expected = expected_ordered_counts[i][j][b];

                g[i][j][b] = (expected > 0.0) ? (observed / expected) : 0.0;
            }
        }
    }

    return g;
}

std::tuple<std::vector<double>, std::vector<std::vector<std::vector<double>>>>
Result::get_uncorrelated_Cij(std::size_t n_bins) const
{
    this->validate_partial_pair_inputs(n_bins);

    const RadialGrid grid = get_radial_grid(n_bins);
    const std::vector<std::size_t> counts = count_particles_per_class();

    const std::size_t K = static_cast<std::size_t>(number_of_classes);
    const std::size_t B = n_bins;

    std::vector<std::vector<std::vector<double>>> expected = make_zero_matrix_expected(K, B);

    const double V = domain->volume;

    for (std::size_t i = 0; i < K; ++i) {
        const double Ni = static_cast<double>(counts[i]);
        if (Ni <= 0.0)
            continue;

        for (std::size_t j = 0; j < K; ++j) {
            const double Nj = static_cast<double>(counts[j]);
            if (Nj <= 0.0)
                continue;

            const double neighbor_count = (i == j) ? std::max(0.0, Nj - 1.0) : Nj;
            const double rho_j = neighbor_count / V;

            for (std::size_t b = 0; b < B; ++b) {
                expected[i][j][b] = Ni * rho_j * grid.shell_volumes[b];
            }
        }
    }

    return { grid.centers, expected };
}

std::tuple<std::vector<double>, std::vector<std::vector<std::vector<double>>>>
Result::compute_partial_pair_correlation_function(
    std::size_t n_bins,
    std::size_t maximum_pairs,
    std::uint64_t random_seed) const
{
    this->validate_partial_pair_inputs(n_bins);

    const std::size_t N = this->sphere_configuration->center_positions.size();

    const RadialGrid grid = get_radial_grid(n_bins);

    std::size_t examined_unordered_pairs = 0;

    const auto histogram_ordered = build_partial_histogram_ordered(
        n_bins,
        grid.r_max,
        grid.dr,
        maximum_pairs,
        random_seed,
        examined_unordered_pairs
    );

    if (examined_unordered_pairs == 0) {
        const std::size_t K = static_cast<std::size_t>(number_of_classes);
        return { grid.centers, make_zero_matrix_double(K, n_bins) };
    }

    const auto uncorrelated = get_uncorrelated_Cij(n_bins);
    const auto& expected_ordered_counts = std::get<1>(uncorrelated);

    const double sampling_scale = compute_pair_sampling_scale(N, examined_unordered_pairs);

    const auto g = compute_g_from_histogram_and_expected(
        histogram_ordered,
        expected_ordered_counts,
        sampling_scale
    );

    return { grid.centers, g };
}
