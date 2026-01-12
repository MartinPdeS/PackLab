#include "analytical/domain_analytical/domain_analytical.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static constexpr double PI = 3.141592653589793238462643383279502884;

static void require(bool condition, const char* message) {
    if (!condition) {
        throw std::invalid_argument(message);
    }
}

PolydisperseDomain::PolydisperseDomain(
    double size_meters_,
    std::vector<double> particle_radii_,
    double volume_fraction_,
    std::vector<double> number_fractions_,
    RoundingMode rounding_mode
)
    : size_meters(size_meters_),
      particle_radii(std::move(particle_radii_)),
      volume_fraction(volume_fraction_),
      number_fractions(std::move(number_fractions_)),
      rounding_mode_(rounding_mode)
{
    require(std::isfinite(size_meters) && size_meters > 0.0, "size must be positive.");
    require(std::isfinite(volume_fraction) && volume_fraction > 0.0 && volume_fraction <= 1.0, "volume_fraction must be in (0, 1].");
    require(!particle_radii.empty(), "particle_radii must not be empty.");

    for (double r : particle_radii) {
        require(std::isfinite(r) && r > 0.0, "particle_radii must be strictly positive and finite.");
    }

    this->validate_and_normalize_number_fractions(number_fractions, particle_radii.size());
}

void PolydisperseDomain::validate_and_normalize_number_fractions(
    std::vector<double>& number_fractions,
    std::size_t expected_size
) {
    require(number_fractions.size() == expected_size, "number_fractions must have the same length as particle_radii.");

    double total = 0.0;
    for (double& v : number_fractions) {
        require(std::isfinite(v), "number_fractions must be finite.");
        require(v >= 0.0, "number_fractions must be non negative.");
        total += v;
    }

    require(total > 0.0, "number_fractions must sum to a positive value.");

    for (double& v : number_fractions) {
        v /= total;
    }
}

double PolydisperseDomain::volume() const {
    return size_meters * size_meters * size_meters;
}

std::vector<double> PolydisperseDomain::get_particle_volumes() const {
    std::vector<double> volumes(particle_radii.size(), 0.0);

    const double prefactor = (4.0 / 3.0) * PI;
    for (std::size_t i = 0; i < particle_radii.size(); ++i) {
        const double r = particle_radii[i];
        volumes[i] = prefactor * r * r * r;
    }
    return volumes;
}

double PolydisperseDomain::get_total_particle_volume() const {
    return volume_fraction * volume();
}

double PolydisperseDomain::get_mean_particle_volume_number_weighted() const {
    const auto volumes = this->get_particle_volumes();
    double accum = 0.0;

    for (std::size_t i = 0; i < volumes.size(); ++i) {
        accum += number_fractions[i] * volumes[i];
    }
    return accum;
}

std::int64_t PolydisperseDomain::apply_rounding_scalar(double value) const {
    if (rounding_mode_ == RoundingMode::Floor) {
        return static_cast<std::int64_t>(std::floor(value));
    }
    return static_cast<std::int64_t>(std::llround(value));
}

std::vector<std::int64_t> PolydisperseDomain::apply_rounding_vector(const std::vector<double>& values) const {
    std::vector<std::int64_t> out(values.size(), 0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        out[i] = this->apply_rounding_scalar(values[i]);
        if (out[i] < 0) {
            out[i] = 0;
        }
    }
    return out;
}

std::int64_t PolydisperseDomain::get_number_of_particles_total() const {
    const double total_occupied_volume = this->get_total_particle_volume();
    const double mean_volume = this->get_mean_particle_volume_number_weighted();

    require(mean_volume > 0.0, "mean particle volume must be positive.");

    const double estimated_total = total_occupied_volume / mean_volume;
    const std::int64_t rounded = this->apply_rounding_scalar(estimated_total);

    return (rounded < 0) ? 0 : rounded;
}

std::vector<std::int64_t> PolydisperseDomain::get_number_of_particles_per_radius() const {
    const std::int64_t total_count = this->get_number_of_particles_total();

    std::vector<double> raw_counts(number_fractions.size(), 0.0);
    for (std::size_t i = 0; i < number_fractions.size(); ++i) {
        raw_counts[i] = number_fractions[i] * static_cast<double>(total_count);
    }

    std::vector<std::int64_t> counts = apply_rounding_vector(raw_counts);

    std::int64_t sum_counts = 0;
    for (std::int64_t c : counts) {
        sum_counts += c;
    }

    std::int64_t difference = total_count - sum_counts;
    if (difference == 0) {
        return counts;
    }

    std::vector<double> residuals(raw_counts.size(), 0.0);
    for (std::size_t i = 0; i < raw_counts.size(); ++i) {
        residuals[i] = raw_counts[i] - static_cast<double>(counts[i]);
    }

    std::vector<std::size_t> indices(residuals.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    if (difference > 0) {
        std::sort(indices.begin(), indices.end(), [&](std::size_t a, std::size_t b) {
            return residuals[a] > residuals[b];
        });

        for (std::int64_t k = 0; k < difference; ++k) {
            counts[indices[static_cast<std::size_t>(k % static_cast<std::int64_t>(indices.size()))]] += 1;
        }
        return counts;
    }

    difference = -difference;

    std::sort(indices.begin(), indices.end(), [&](std::size_t a, std::size_t b) {
        return residuals[a] < residuals[b];
    });

    std::int64_t removed = 0;
    for (std::size_t idx : indices) {
        if (removed >= difference) {
            break;
        }
        if (counts[idx] > 0) {
            counts[idx] -= 1;
            removed += 1;
        }
    }

    return counts;
}

std::vector<double> PolydisperseDomain::get_particle_densities_per_radius() const {
    const double V = volume();
    const auto counts = this->get_number_of_particles_per_radius();

    std::vector<double> densities(counts.size(), 0.0);
    for (std::size_t i = 0; i < counts.size(); ++i) {
        densities[i] = static_cast<double>(counts[i]) / V;
    }
    return densities;
}

double PolydisperseDomain::get_particle_density_total() const {
    return static_cast<double>(this->get_number_of_particles_total()) / volume();
}

std::vector<double> PolydisperseDomain::get_volume_fraction_per_radius() const {
    const double V = volume();
    const auto counts = this->get_number_of_particles_per_radius();
    const auto volumes = this->get_particle_volumes();

    std::vector<double> vf(counts.size(), 0.0);
    for (std::size_t i = 0; i < counts.size(); ++i) {
        vf[i] = (static_cast<double>(counts[i]) * volumes[i]) / V;
    }
    return vf;
}

std::vector<double> PolydisperseDomain::sample_particle_radii(
    std::int64_t number_of_samples,
    std::uint64_t seed
) const {
    require(number_of_samples > 0, "number_of_samples must be > 0.");

    std::mt19937_64 rng(seed);

    std::discrete_distribution<std::size_t> pick(number_fractions.begin(), number_fractions.end());

    std::vector<double> out(static_cast<std::size_t>(number_of_samples), 0.0);
    for (std::int64_t i = 0; i < number_of_samples; ++i) {
        const std::size_t index = pick(rng);
        out[static_cast<std::size_t>(i)] = particle_radii[index];
    }
    return out;
}


static std::string format_g(double value, int precision) {
    std::ostringstream oss;
    oss.setf(std::ios::fmtflags(0), std::ios::floatfield);
    oss << std::setprecision(precision) << std::defaultfloat << value;
    return oss.str();
}

static std::string pad_right(const std::string& s, std::size_t width) {
    if (s.size() >= width) {
        return s;
    }
    return s + std::string(width - s.size(), ' ');
}

static std::string pad_left(const std::string& s, std::size_t width) {
    if (s.size() >= width) {
        return s;
    }
    return std::string(width - s.size(), ' ') + s;
}

std::string PolydisperseDomain::bins_table(int precision) const {


    const std::size_t n_bins = this->particle_radii.size();


    const std::vector<double> volumes_m3 = this->get_particle_volumes();
    const std::vector<std::int64_t> counts_i64 = this->get_number_of_particles_per_radius();
    const std::vector<double> densities = this->get_particle_densities_per_radius();
    const std::vector<double> vf = this->get_volume_fraction_per_radius();


    std::vector<std::string> headers = {
        "bin",
        "radius (m)",
        "particle_volume (m^3)",
        "number_fraction",
        "volume_fraction",
        "particle_count",
        "number_density (1/m^3)",
    };

    std::vector<std::vector<std::string>> rows;
    rows.reserve(n_bins);

    for (std::size_t i = 0; i < n_bins; ++i) {
        std::vector<std::string> row;
        row.reserve(headers.size());

        row.push_back(std::to_string(static_cast<int>(i)));
        row.push_back(format_g(this->particle_radii[i], precision));
        row.push_back(format_g(volumes_m3[i], precision));
        row.push_back(format_g(this->number_fractions[i], precision));
        row.push_back(format_g(vf[i], precision));
        row.push_back(std::to_string(static_cast<long long>(counts_i64[i])));
        row.push_back(format_g(densities[i], precision));

        rows.push_back(std::move(row));
    }


    std::vector<std::size_t> widths(headers.size(), 0);
    for (std::size_t c = 0; c < headers.size(); ++c) {
        widths[c] = headers[c].size();
    }
    for (const auto& row : rows) {
        for (std::size_t c = 0; c < row.size(); ++c) {
            widths[c] = std::max(widths[c], row[c].size());
        }
    }

    auto make_separator = [&]() {
        std::ostringstream oss;
        oss << "+";
        for (std::size_t c = 0; c < widths.size(); ++c) {
            oss << std::string(widths[c] + 2, '-') << "+";
        }
        return oss.str();
    };

    const std::string sep = make_separator();

    std::ostringstream out;

    out << sep << "\n";
    out << "|";
    for (std::size_t c = 0; c < headers.size(); ++c) {
        out << " " << pad_right(headers[c], widths[c]) << " |";
    }
    out << "\n";
    out << sep << "\n";

    for (const auto& row : rows) {
        out << "|";
        for (std::size_t c = 0; c < row.size(); ++c) {
            const bool right_align =
                (c == 0) || (c == 5) || (c == 1) || (c == 2) || (c == 3) || (c == 4) || (c == 6);

            if (right_align) {
                out << " " << pad_left(row[c], widths[c]) << " |";
            } else {
                out << " " << pad_right(row[c], widths[c]) << " |";
            }
        }
        out << "\n";
    }

    out << sep << "\n";
    return out.str();
}

void PolydisperseDomain::print_bins(int precision) const {
    std::cout << bins_table(precision);
}
