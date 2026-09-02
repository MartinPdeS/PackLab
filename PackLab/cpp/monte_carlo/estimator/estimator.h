#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "radius_sampler/discrete.h"
#include "radius_sampler/normal.h"
#include "radius_sampler/log_normal.h"
#include "radius_sampler/uniform.h"
#include "radius_sampler/constant.h"
#include "monte_carlo/domain/domain.h"
#include "monte_carlo/simulator/simulator.h"
#include "monte_carlo/utils/utils.h"
#include "monte_carlo/utils/numpy.h"

struct EstimateResult {
    std::shared_ptr<MCDomain> domain;
    std::vector<double> centers;   // size B
    std::vector<std::vector<std::vector<double>>> mean_g; // size K x K x B
    std::vector<std::vector<std::vector<double>>> std_g;  // size K x K x B
    std::size_t number_of_species = 0; // K
    std::size_t number_of_bins = 0;    // B
};

struct EstimatorStatistics {
    std::size_t requested_samples = 0;
    std::size_t completed_samples = 0;
    std::size_t attempted_insertions = 0;
    std::size_t accepted_insertions = 0;
    std::size_t rejected_insertions = 0;
    std::size_t total_spheres = 0;
    double total_runtime_seconds = 0.0;
    double mean_packing_fraction = 0.0;

    double acceptance_rate() const {
        if (attempted_insertions == 0) return 0.0;
        return static_cast<double>(accepted_insertions) /
            static_cast<double>(attempted_insertions);
    }

    double mean_sphere_count() const {
        if (completed_samples == 0) return 0.0;
        return static_cast<double>(total_spheres) /
            static_cast<double>(completed_samples);
    }

    double mean_runtime_seconds() const {
        if (completed_samples == 0) return 0.0;
        return total_runtime_seconds / static_cast<double>(completed_samples);
    }

    void print() const {
        const auto print_row = [](const std::string& key, const std::string& value) {
            std::cout << "| " << std::left << std::setw(31) << key
                      << " | " << std::right << std::setw(16) << value << " |\n";
        };

        std::cout << "+---------------------------------+------------------+\n";
        print_row("completed samples", std::to_string(completed_samples) + "/" + std::to_string(requested_samples));
        print_row("attempted insertions", std::to_string(attempted_insertions));
        print_row("accepted insertions", std::to_string(accepted_insertions));
        print_row("rejected insertions", std::to_string(rejected_insertions));
        print_row("acceptance rate", std::to_string(acceptance_rate()));
        print_row("mean sphere count", std::to_string(mean_sphere_count()));
        print_row("mean packing fraction", std::to_string(mean_packing_fraction));
        print_row("total runtime [seconds]", std::to_string(total_runtime_seconds));
        print_row("mean runtime [seconds]", std::to_string(mean_runtime_seconds()));
        std::cout << "+---------------------------------+------------------+\n" << std::endl;
    }
};

class Estimator {
public:
    Estimator(
        const std::shared_ptr<MCDomain>& domain,
        const std::shared_ptr<RadiusSampler>& radius_sampler,
        const std::shared_ptr<Options>& options,
        std::size_t number_of_bins
    )
        : domain(domain),
          radius_sampler(radius_sampler),
          options(options),
          number_of_bins(number_of_bins) {}

    EstimateResult estimate(
        const std::size_t number_of_samples,
        const std::size_t maximum_pairs = 0,
        const bool progress = false,
        const std::size_t progress_interval = 1
    ) {
        if (number_of_samples == 0) {
            throw std::invalid_argument("number_of_samples must be > 0");
        }
        if (progress_interval == 0) {
            throw std::invalid_argument("progress_interval must be > 0.");
        }

        statistics = EstimatorStatistics{};
        statistics.requested_samples = number_of_samples;

        if (progress) {
            std::cout << "PackingEstimator progress\n"
                      << "  " << std::right << std::setw(10) << "sample"
                      << "  " << std::setw(10) << "accepted"
                      << "  " << std::setw(10) << "attempted"
                      << "  " << std::setw(15) << "acceptance rate"
                      << "  " << std::setw(16) << "packing fraction"
                      << std::endl;
        }

        Simulator simulator(domain, radius_sampler, options);

        std::vector<double> centers;
        std::vector<std::vector<std::vector<double>>> mean_g;
        std::vector<std::vector<std::vector<double>>> m2_g;

        std::size_t sample_count = 0;
        std::size_t inferred_number_of_bins = 0;
        std::size_t inferred_number_of_species = 0;

        for (std::size_t sample_index = 0; sample_index < number_of_samples; ++sample_index) {
            simulator.reset();
            Result result = simulator.run();

            const Statistics& sample_statistics = result.statistics;
            statistics.completed_samples += 1;
            statistics.attempted_insertions += sample_statistics.attempted_insertions;
            statistics.accepted_insertions += sample_statistics.accepted_insertions;
            statistics.rejected_insertions += sample_statistics.rejected_insertions;
            statistics.total_spheres += sample_statistics.sphere_count;
            statistics.total_runtime_seconds += sample_statistics.total_runtime_seconds;
            statistics.mean_packing_fraction +=
                (sample_statistics.packing_fraction_geometry - statistics.mean_packing_fraction) /
                static_cast<double>(statistics.completed_samples);
            if (progress && (
                statistics.completed_samples == 1 ||
                statistics.completed_samples % progress_interval == 0 ||
                statistics.completed_samples == number_of_samples
            )) {
                const double acceptance_rate = sample_statistics.attempted_insertions == 0
                    ? 0.0
                    : static_cast<double>(sample_statistics.accepted_insertions) /
                        static_cast<double>(sample_statistics.attempted_insertions);

                std::cout << "  " << std::right << std::setw(10)
                          << std::to_string(statistics.completed_samples) + "/" + std::to_string(number_of_samples)
                          << "  " << std::setw(10) << sample_statistics.accepted_insertions
                          << "  " << std::setw(10) << sample_statistics.attempted_insertions
                          << "  " << std::setw(14) << std::fixed << std::setprecision(3)
                          << 100.0 * acceptance_rate << "%"
                          << "  " << std::setw(16) << std::setprecision(6)
                          << sample_statistics.packing_fraction_geometry
                          << std::endl;
            }

            auto [centers_i, g_matrix_i] =
                result.compute_partial_pair_correlation_function(number_of_bins, maximum_pairs);

            sample_count += 1;

            if (sample_count == 1) {
                centers = std::move(centers_i);
                inferred_number_of_bins = centers.size();

                if (inferred_number_of_bins == 0) {
                    throw std::runtime_error("compute_partial_pair_correlation_function returned empty centers");
                }

                inferred_number_of_species = g_matrix_i.size();
                if (inferred_number_of_species == 0) {
                    throw std::runtime_error("compute_partial_pair_correlation_function returned empty g_matrix");
                }

                for (std::size_t i = 0; i < inferred_number_of_species; ++i) {
                    if (g_matrix_i[i].size() != inferred_number_of_species) {
                        throw std::runtime_error("g_matrix is not K x K along the first two dimensions");
                    }
                    for (std::size_t j = 0; j < inferred_number_of_species; ++j) {
                        if (g_matrix_i[i][j].size() != inferred_number_of_bins) {
                            throw std::runtime_error("g_matrix bin dimension does not match centers size");
                        }
                    }
                }

                mean_g = std::move(g_matrix_i);
                m2_g = make_zero_accumulator_like(mean_g);
                continue;
            }

            if (centers_i.size() != inferred_number_of_bins) {
                throw std::runtime_error("centers size changed across samples");
            }
            if (g_matrix_i.size() != inferred_number_of_species) {
                throw std::runtime_error("g_matrix K dimension changed across samples");
            }

            for (std::size_t i = 0; i < inferred_number_of_species; ++i) {
                if (g_matrix_i[i].size() != inferred_number_of_species) {
                    throw std::runtime_error("g_matrix second dimension changed across samples");
                }
                for (std::size_t j = 0; j < inferred_number_of_species; ++j) {
                    if (g_matrix_i[i][j].size() != inferred_number_of_bins) {
                        throw std::runtime_error("g_matrix bin dimension changed across samples");
                    }
                }
            }

            const double inv_count = 1.0 / static_cast<double>(sample_count);

            // Welford update for every element: (i, j, b)
            for (std::size_t i = 0; i < inferred_number_of_species; ++i) {
                for (std::size_t j = 0; j < inferred_number_of_species; ++j) {
                    for (std::size_t b = 0; b < inferred_number_of_bins; ++b) {
                        const double x = g_matrix_i[i][j][b];
                        const double delta = x - mean_g[i][j][b];
                        mean_g[i][j][b] += delta * inv_count;
                        const double delta2 = x - mean_g[i][j][b];
                        m2_g[i][j][b] += delta * delta2;
                    }
                }
            }
        }

        std::vector<std::vector<std::vector<double>>> std_g = make_zero_accumulator_like(mean_g);

        if (number_of_samples > 1) {
            const double inv = 1.0 / static_cast<double>(number_of_samples - 1);
            for (std::size_t i = 0; i < inferred_number_of_species; ++i) {
                for (std::size_t j = 0; j < inferred_number_of_species; ++j) {
                    for (std::size_t b = 0; b < inferred_number_of_bins; ++b) {
                        const double variance = m2_g[i][j][b] * inv;
                        std_g[i][j][b] = std::sqrt(variance);
                    }
                }
            }
        }

        EstimateResult out;
        out.centers = std::move(centers);

        out.domain = domain;

        out.mean_g = std::move(mean_g);
        out.std_g = std::move(std_g);
        out.number_of_species = inferred_number_of_species;
        out.number_of_bins = inferred_number_of_bins;
        return out;
    }

private:
    static std::vector<std::vector<std::vector<double>>>
    make_zero_accumulator_like(const std::vector<std::vector<std::vector<double>>>& reference) {
        std::vector<std::vector<std::vector<double>>> accumulator(reference.size());
        for (std::size_t i = 0; i < reference.size(); ++i) {
            accumulator[i].resize(reference[i].size());
            for (std::size_t j = 0; j < reference[i].size(); ++j) {
                accumulator[i][j].assign(reference[i][j].size(), 0.0);
            }
        }
        return accumulator;
    }

private:
    std::shared_ptr<MCDomain> domain;
    std::shared_ptr<RadiusSampler> radius_sampler;
    std::shared_ptr<Options> options;
    std::size_t number_of_bins = 0;

public:
    EstimatorStatistics statistics;
};
