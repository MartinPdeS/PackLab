#include <cstddef>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "monte_carlo/domain/domain.h"
#include "monte_carlo/radius_sampler/radius_sampler.h"
#include "monte_carlo/simulator/simulator.h"
#include "monte_carlo/utils/utils.h"
#include "monte_carlo/utils/numpy.h"

struct EstimateResult {
    std::shared_ptr<Domain> domain;
    std::vector<double> centers;   // size B
    std::vector<std::vector<std::vector<double>>> mean_g; // size K x K x B
    std::vector<std::vector<std::vector<double>>> std_g;  // size K x K x B
    std::size_t number_of_species = 0; // K
    std::size_t number_of_bins = 0;    // B
};

class Estimator {
public:
    Estimator(
        const std::shared_ptr<Domain>& domain,
        const std::shared_ptr<RadiusSampler>& radius_sampler,
        const std::shared_ptr<Options>& options,
        std::size_t number_of_bins
    )
        : domain(domain),
          radius_sampler(radius_sampler),
          options(options),
          number_of_bins(number_of_bins) {}

    EstimateResult estimate(const std::size_t number_of_samples, const std::size_t maximum_pairs = 0) const {
        if (number_of_samples == 0) {
            throw std::invalid_argument("number_of_samples must be > 0");
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
    std::shared_ptr<Domain> domain;
    std::shared_ptr<RadiusSampler> radius_sampler;
    std::shared_ptr<Options> options;
    std::size_t number_of_bins = 0;
};
