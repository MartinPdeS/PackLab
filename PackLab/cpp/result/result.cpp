#include "result.h"

namespace {

// Parameters used internally during g(r) computation
struct PairCorrelationParameters {
    double maximum_distance = 0.0;
    double radial_bin_width = 0.0;
    bool periodic = true;
};

// Resolve grid cell size for the grid-based histogram
double resolve_grid_cell_size(
    double grid_cell_size,
    double maximum_distance,
    const Domain& domain
) {
    const double Lx = domain.length_x();
    const double Ly = domain.length_y();
    const double Lz = domain.length_z();
    const bool periodic = domain.use_periodic_boundaries();

    const double periodic_limit = 0.5 * std::min({Lx, Ly, Lz});
    const double nonperiodic_limit =
        std::sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

    const double max_allowed =
        periodic ? periodic_limit : nonperiodic_limit;

    double resolved = grid_cell_size;

    if (resolved <= 0.0) {
        resolved = maximum_distance;
    } else {
        const double min_reasonable = maximum_distance / 4.0;
        if (resolved < min_reasonable) resolved = min_reasonable;
        if (resolved > max_allowed)    resolved = max_allowed;
    }

    return resolved;
}

} // namespace



// -----------------------------------------------------------------------------
// Monte Carlo histogram fill
// -----------------------------------------------------------------------------
void Result::fill_pair_correlation_histogram_monte_carlo(
    std::vector<std::size_t>& histogram,
    std::size_t bins,
    std::size_t maximum_number_of_pairs,
    std::uint64_t random_seed,
    double maximum_distance,
    double radial_bin_width
) const {
    const std::size_t N = this->sphere_configuration.center_positions_values_.size();
    histogram.assign(bins, 0);

    std::mt19937_64 rng(random_seed);

    for (std::size_t k = 0; k < maximum_number_of_pairs; ++k) {
        const double r = sample_distance(
            N,
            this->sphere_configuration.center_positions_values_,
            rng,
            domain
        );

        if (r >= maximum_distance) continue;

        const std::size_t index = static_cast<std::size_t>(r / radial_bin_width);
        if (index < bins) histogram[index] += 1;
    }
}



// -----------------------------------------------------------------------------
// Spatial grid histogram fill
// -----------------------------------------------------------------------------
void Result::fill_pair_correlation_histogram_grid(
    std::vector<std::size_t>& histogram,
    std::size_t bins,
    double maximum_distance,
    double radial_bin_width,
    double grid_cell_size
) const {
    const std::size_t N = this->sphere_configuration.center_positions_values_.size();
    histogram.assign(bins, 0);

    const double resolved_cell_size =
        resolve_grid_cell_size(grid_cell_size, maximum_distance, domain);

    SpatialGridIndex grid(resolved_cell_size, domain);

    // Populate grid
    for (std::size_t i = 0; i < N; ++i)
        grid.insert_sphere(i, this->sphere_configuration.center_positions_values_[i]);

    const std::int64_t cell_range =
        static_cast<std::int64_t>(std::ceil(maximum_distance / grid.cell_size()));

    const bool periodic = domain.use_periodic_boundaries();
    const double Lx = domain.length_x();
    const double Ly = domain.length_y();
    const double Lz = domain.length_z();

    std::vector<std::size_t> neighbors;
    neighbors.reserve(512);

    for (std::size_t i = 0; i < N; ++i) {
        neighbors.clear();
        grid.append_neighbor_sphere_indices(
            this->sphere_configuration.center_positions_values_[i],
            cell_range,
            neighbors
        );

        const Vector3d& pi = this->sphere_configuration.center_positions_values_[i];

        for (std::size_t j : neighbors) {
            if (j <= i) continue;

            const Vector3d& pj = this->sphere_configuration.center_positions_values_[j];

            double dx = pj.x - pi.x;
            double dy = pj.y - pi.y;
            double dz = pj.z - pi.z;

            if (periodic) {
                dx = domain.minimum_image_displacement(dx, Lx);
                dy = domain.minimum_image_displacement(dy, Ly);
                dz = domain.minimum_image_displacement(dz, Lz);
            }

            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r >= maximum_distance) continue;

            const std::size_t index = static_cast<std::size_t>(r / radial_bin_width);
            if (index < bins) histogram[index] += 1;
        }
    }
}



// -----------------------------------------------------------------------------
// Minimum image helper
// -----------------------------------------------------------------------------
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



// -----------------------------------------------------------------------------
// Random pair distance sampler
// -----------------------------------------------------------------------------
double Result::sample_distance(
    std::size_t number_of_particles,
    const std::vector<Vector3d>& positions,
    std::mt19937_64& rng,
    const Domain& domain
) const {
    std::uniform_int_distribution<std::size_t> pick(0, number_of_particles - 1);

    std::size_t i = 0, j = 0;
    while ((i = pick(rng)) == (j = pick(rng))) {}

    Vector3d d = positions[j] - positions[i];

    if (domain.use_periodic_boundaries())
        d = apply_minimum_image(d, domain);

    return d.norm();
}



// -----------------------------------------------------------------------------
// Top-level dispatcher for g(r)
// -----------------------------------------------------------------------------
void Result::compute_pair_correlation_function(
    std::size_t bins,
    double maximum_distance,
    const std::string& method,
    std::size_t maximum_number_of_pairs,
    std::uint64_t random_seed,
    double grid_cell_size
) {
    pair_correlation_centers_.clear();
    pair_correlation_values_.clear();

    const std::size_t N = this->sphere_configuration.center_positions_values_.size();
    if (N < 2 || bins == 0)
        return;

    const double Lx = domain.length_x();
    const double Ly = domain.length_y();
    const double Lz = domain.length_z();
    const bool periodic = domain.use_periodic_boundaries();

    const double periodic_limit = 0.5 * std::min({Lx, Ly, Lz});
    const double nonperiodic_limit =
        std::sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

    const double default_max_dist =
        periodic ? periodic_limit : nonperiodic_limit;

    double r_max = maximum_distance <= 0.0
        ? default_max_dist
        : maximum_distance;

    if (periodic && r_max > periodic_limit) r_max = periodic_limit;
    if (!periodic && r_max > nonperiodic_limit) r_max = nonperiodic_limit;

    const double dr = r_max / static_cast<double>(bins);

    pair_correlation_centers_.resize(bins);
    for (std::size_t i = 0; i < bins; ++i) {
        pair_correlation_centers_[i] =
            (static_cast<double>(i) + 0.5) * dr;
    }

    // Normalize method name
    std::string m = method;
    std::transform(m.begin(), m.end(), m.begin(),
                   [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });

    if (m == "grid") {
        pair_correlation_values_ =
            compute_pair_correlation_values_once_grid(
                bins, r_max, grid_cell_size
            );
        return;
    }

    if (m == "monte-carlo" || m == "mc") {
        if (maximum_number_of_pairs == 0)
            throw std::invalid_argument(
                "maximum_number_of_pairs must be nonzero for Monte Carlo."
            );

        pair_correlation_values_ =
            compute_pair_correlation_values_once(
                bins, maximum_number_of_pairs, random_seed, r_max
            );
        return;
    }

    throw std::invalid_argument("method must be 'grid' or 'monte-carlo'.");
}



// -----------------------------------------------------------------------------
// Monte Carlo g(r) computation
// -----------------------------------------------------------------------------
std::vector<double> Result::compute_pair_correlation_values_once(
    std::size_t bins,
    std::size_t maximum_number_of_pairs,
    std::uint64_t random_seed,
    double maximum_distance
) const {
    const std::size_t N = this->sphere_configuration.center_positions_values_.size();
    std::vector<double> values(bins, 0.0);

    if (N < 2 || bins == 0 || maximum_number_of_pairs == 0)
        return values;

    const double V = domain.volume();
    const double dr = maximum_distance / static_cast<double>(bins);

    std::vector<std::size_t> histogram(bins, 0);

    std::mt19937_64 rng(random_seed);
    std::uniform_int_distribution<std::size_t> pick(0, N - 1);

    const bool periodic = domain.use_periodic_boundaries();
    const double Lx = domain.length_x();
    const double Ly = domain.length_y();
    const double Lz = domain.length_z();

    for (std::size_t k = 0; k < maximum_number_of_pairs; ++k) {
        std::size_t i = pick(rng);
        std::size_t j = pick(rng);
        while (j == i) j = pick(rng);

        const Vector3d& pi = this->sphere_configuration.center_positions_values_[i];
        const Vector3d& pj = this->sphere_configuration.center_positions_values_[j];

        double dx = pj.x - pi.x;
        double dy = pj.y - pi.y;
        double dz = pj.z - pi.z;

        if (periodic) {
            dx = domain.minimum_image_displacement(dx, Lx);
            dy = domain.minimum_image_displacement(dy, Ly);
            dz = domain.minimum_image_displacement(dz, Lz);
        }

        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (r >= maximum_distance) continue;

        const std::size_t b = static_cast<std::size_t>(r / dr);
        if (b < bins) histogram[b] += 1;
    }

    const double M = static_cast<double>(maximum_number_of_pairs);

    for (std::size_t b = 0; b < bins; ++b) {
        const double r0 = static_cast<double>(b) * dr;
        const double r1 = static_cast<double>(b + 1) * dr;

        const double shell_vol =
            (4.0 * PI / 3.0) * (std::pow(r1, 3) - std::pow(r0, 3));

        const double expected = M * (shell_vol / V);
        values[b] = expected > 0.0
            ? static_cast<double>(histogram[b]) / expected
            : 0.0;
    }

    return values;
}



// -----------------------------------------------------------------------------
// Grid-based g(r)
// -----------------------------------------------------------------------------
std::vector<double> Result::compute_pair_correlation_values_once_grid(
    std::size_t bins,
    double maximum_distance,
    double grid_cell_size
) const {
    const std::size_t N = this->sphere_configuration.center_positions_values_.size();
    std::vector<double> values(bins, 0.0);

    if (N < 2 || bins == 0)
        return values;

    const double V = domain.volume();

    std::vector<std::size_t> histogram(bins, 0);

    const double Lx = domain.length_x();
    const double Ly = domain.length_y();
    const double Lz = domain.length_z();
    const bool periodic = domain.use_periodic_boundaries();

    const double periodic_limit = 0.5 * std::min({Lx, Ly, Lz});
    const double nonperiodic_limit = std::sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    const double default_max_dist = periodic ? periodic_limit : nonperiodic_limit;

    double r_max = maximum_distance <= 0.0
        ? default_max_dist
        : maximum_distance;

    if (periodic && r_max > periodic_limit) r_max = periodic_limit;
    if (!periodic && r_max > nonperiodic_limit) r_max = nonperiodic_limit;

    const double dr = r_max / static_cast<double>(bins);

    SpatialGridIndex grid(
        resolve_grid_cell_size(grid_cell_size, r_max, domain),
        domain
    );

    for (std::size_t i = 0; i < N; ++i)
        grid.insert_sphere(i, this->sphere_configuration.center_positions_values_[i]);

    const std::int64_t cell_range =
        static_cast<std::int64_t>(std::ceil(r_max / grid.cell_size()));

    std::vector<std::size_t> neighbors;
    neighbors.reserve(512);

    for (std::size_t i = 0; i < N; ++i) {
        neighbors.clear();
        grid.append_neighbor_sphere_indices(
            this->sphere_configuration.center_positions_values_[i],
            cell_range,
            neighbors
        );

        const Vector3d& pi = this->sphere_configuration.center_positions_values_[i];

        for (std::size_t j : neighbors) {
            if (j <= i) continue;

            const Vector3d& pj = this->sphere_configuration.center_positions_values_[j];

            double dx = pj.x - pi.x;
            double dy = pj.y - pi.y;
            double dz = pj.z - pi.z;

            if (periodic) {
                dx = domain.minimum_image_displacement(dx, Lx);
                dy = domain.minimum_image_displacement(dy, Ly);
                dz = domain.minimum_image_displacement(dz, Lz);
            }

            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r >= r_max) continue;

            const std::size_t b = static_cast<std::size_t>(r / dr);
            if (b < bins) histogram[b] += 1;
        }
    }

    const double M =
        static_cast<double>(N) * static_cast<double>(N - 1) * 0.5;

    for (std::size_t b = 0; b < bins; ++b) {
        const double r0 = static_cast<double>(b) * dr;
        const double r1 = static_cast<double>(b + 1) * dr;

        const double shell_vol =
            (4.0 * PI / 3.0) * (std::pow(r1, 3) - std::pow(r0, 3));

        const double expected = M * (shell_vol / V);
        values[b] = expected > 0.0
            ? static_cast<double>(histogram[b]) / expected
            : 0.0;
    }

    return values;
}



// -----------------------------------------------------------------------------
// Partial g_ij(r) for class-labeled polydisperse system
// -----------------------------------------------------------------------------
std::tuple<
    std::vector<double>,
    std::vector<std::vector<std::vector<double>>>
>
Result::compute_partial_pair_correlation_function(
    std::size_t number_of_distance_bins,
    std::size_t maximum_pairs,
    std::uint64_t random_seed
) const {
    const std::size_t N = this->sphere_configuration.center_positions_values_.size();
    if (N < 2)
        return { {}, {} };

    if (this->sphere_configuration.class_index_values_.empty())
        throw std::runtime_error(
            "Partial pair correlation requires this->sphere_configuration.class_index_values_ to be set."
        );

    const int K = number_of_classes_;
    if (K <= 0)
        throw std::runtime_error(
            "number_of_classes_ must be > 0."
        );

    // Count population of each class
    std::vector<std::size_t> count_per_class(K, 0);
    for (int c : this->sphere_configuration.class_index_values_) {
        if (c < 0 || c >= K)
            throw std::runtime_error("Invalid class index in this->sphere_configuration.class_index_values_.");
        count_per_class[c] += 1;
    }

    // Domain geometry
    const double Lx = domain.length_x();
    const double Ly = domain.length_y();
    const double Lz = domain.length_z();
    const double V  = domain.volume();

    // Density of each class
    std::vector<double> rho(K);
    for (int j = 0; j < K; ++j)
        rho[j] = static_cast<double>(count_per_class[j]) / V;

    // Maximum possible distance
    const double r_max =
        0.5 * std::min({Lx, Ly, Lz});

    const double dr = r_max / static_cast<double>(number_of_distance_bins);

    // histogram[i][j][b]
    std::vector<std::vector<std::vector<std::size_t>>> histogram(
        K,
        std::vector<std::vector<std::size_t>>(
            K,
            std::vector<std::size_t>(number_of_distance_bins, 0)
        )
    );

    std::mt19937_64 rng(random_seed);
    std::uniform_int_distribution<std::size_t> pick(0, N - 1);

    std::size_t sampled = 0;

    const bool periodic = domain.use_periodic_boundaries();

    while (sampled < maximum_pairs) {
        std::size_t p = pick(rng);
        std::size_t q = pick(rng);
        if (p == q) continue;

        int ci = this->sphere_configuration.class_index_values_[p];
        int cj = this->sphere_configuration.class_index_values_[q];

        Vector3d d = this->sphere_configuration.center_positions_values_[q] - this->sphere_configuration.center_positions_values_[p];

        if (periodic)
            d = apply_minimum_image(d, domain);

        const double r = d.norm();
        if (r >= r_max) continue;

        const std::size_t b = static_cast<std::size_t>(r / dr);
        if (b < number_of_distance_bins) {
            histogram[ci][cj][b] += 1;
            sampled += 1;
        }
    }

    // Prepare output
    std::vector<double> centers(number_of_distance_bins);
    for (std::size_t b = 0; b < number_of_distance_bins; ++b)
        centers[b] = (b + 0.5) * dr;

    // g_ij[r][i][j]
    std::vector<std::vector<std::vector<double>>> g_matrix(
        K,
        std::vector<std::vector<double>>(
            K,
            std::vector<double>(number_of_distance_bins, 0.0)
        )
    );

    // Convert histogram -> g_ij(r)
    for (int i = 0; i < K; ++i) {
        const double Ni = static_cast<double>(count_per_class[i]);
        if (Ni <= 0.0) continue;

        for (int j = 0; j < K; ++j) {
            const double rho_j = rho[j];
            if (rho_j <= 0.0) continue;

            for (std::size_t b = 0; b < number_of_distance_bins; ++b) {
                const double r0 = static_cast<double>(b) * dr;
                const double r1 = static_cast<double>(b + 1) * dr;

                const double shell_vol =
                    (4.0 * PI / 3.0) *
                    (std::pow(r1, 3) - std::pow(r0, 3));

                const double expected =
                    Ni * rho_j * shell_vol;

                g_matrix[i][j][b] =
                    expected > 0.0
                    ? static_cast<double>(histogram[i][j][b]) / expected
                    : 0.0;
            }
        }
    }

    return { centers, g_matrix };
}
