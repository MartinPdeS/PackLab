#include "py_solver.h"

#include <sstream>

static constexpr double pi_value = 3.141592653589793238462643383279502884;

static void require(bool condition, const char* message) {
    if (!condition) {
        throw std::invalid_argument(message);
    }
}

static inline double safe_sin_over_x(double x, double small = 1e-12) {
    const double ax = std::abs(x);
    if (ax <= small) {
        return 1.0;
    }
    return std::sin(x) / x;
}

static inline double safe_sin_over_x3_minus_cos_over_x2(double x, double small = 1e-12) {
    const double ax = std::abs(x);
    if (ax <= small) {
        return 1.0 / 3.0;
    }
    return std::sin(x) / (x * x * x) - std::cos(x) / (x * x);
}

// Row major indexing helpers
static inline std::size_t index_2d(std::size_t i, std::size_t j, std::size_t ncols) {
    return i * ncols + j;
}

static inline std::size_t index_3d(std::size_t i, std::size_t j, std::size_t k, std::size_t n_j, std::size_t n_k) {
    return (i * n_j + j) * n_k + k;
}

// Simple LU with partial pivoting, solves A X = B for X, where B has N columns
static void lu_solve_in_place_multiple_rhs(
    std::vector<double>& A,            // size N*N, overwritten by LU
    std::vector<double>& B,            // size N*N, overwritten by solution X
    std::size_t N
) {
    std::vector<std::size_t> piv(N);
    for (std::size_t i = 0; i < N; ++i) {
        piv[i] = i;
    }

    for (std::size_t k = 0; k < N; ++k) {
        // Pivot
        std::size_t pivot_row = k;
        double max_abs = std::abs(A[index_2d(k, k, N)]);
        for (std::size_t i = k + 1; i < N; ++i) {
            const double v = std::abs(A[index_2d(i, k, N)]);
            if (v > max_abs) {
                max_abs = v;
                pivot_row = i;
            }
        }
        require(max_abs > 0.0, "Singular matrix in LU solve.");

        if (pivot_row != k) {
            for (std::size_t j = 0; j < N; ++j) {
                std::swap(A[index_2d(k, j, N)], A[index_2d(pivot_row, j, N)]);
                std::swap(B[index_2d(k, j, N)], B[index_2d(pivot_row, j, N)]);
            }
            std::swap(piv[k], piv[pivot_row]);
        }

        const double Akk = A[index_2d(k, k, N)];
        for (std::size_t i = k + 1; i < N; ++i) {
            A[index_2d(i, k, N)] /= Akk;
            const double Lik = A[index_2d(i, k, N)];
            for (std::size_t j = k + 1; j < N; ++j) {
                A[index_2d(i, j, N)] -= Lik * A[index_2d(k, j, N)];
            }
        }
    }

    // Forward substitution: solve L Y = B
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            double sum = B[index_2d(i, j, N)];
            for (std::size_t k = 0; k < i; ++k) {
                sum -= A[index_2d(i, k, N)] * B[index_2d(k, j, N)];
            }
            B[index_2d(i, j, N)] = sum;
        }
    }

    // Back substitution: solve U X = Y
    for (std::size_t i = N; i-- > 0;) {
        const double Uii = A[index_2d(i, i, N)];
        for (std::size_t j = 0; j < N; ++j) {
            double sum = B[index_2d(i, j, N)];
            for (std::size_t k = i + 1; k < N; ++k) {
                sum -= A[index_2d(i, k, N)] * B[index_2d(k, j, N)];
            }
            B[index_2d(i, j, N)] = sum / Uii;
        }
    }
}

PercusYevickSolver::PercusYevickSolver(
    std::vector<double> densities_per_m3,
    std::vector<double> radii_m,
    std::vector<double> wavenumber_per_m
)
    : densities_per_m3_(std::move(densities_per_m3)),
      radii_m_(std::move(radii_m)),
      wavenumber_per_m_(std::move(wavenumber_per_m))
{
    require(!densities_per_m3_.empty(), "densities must not be empty.");
    require(densities_per_m3_.size() == radii_m_.size(), "densities and radii must have the same length.");
    require(!wavenumber_per_m_.empty(), "wavenumber must not be empty.");

    for (double r : radii_m_) {
        require(std::isfinite(r) && r > 0.0, "radii must be positive and finite.");
    }
    for (double n : densities_per_m3_) {
        require(std::isfinite(n) && n >= 0.0, "densities must be non negative and finite.");
    }
    for (double wavenumber : wavenumber_per_m_) {
        require(std::isfinite(wavenumber) && wavenumber >= 0.0, "wavenumber must be non-negative and finite.");
    }
    for (std::size_t index = 1; index < wavenumber_per_m_.size(); ++index) {
        require(wavenumber_per_m_[index] > wavenumber_per_m_[index - 1], "wavenumber must be strictly increasing.");
    }
}

PercusYevickSolver::PercusYevickSolver(
    std::vector<double> densities_per_m3,
    std::vector<double> radii_m,
    double automatic_radial_resolution_m,
    int automatic_samples_per_oscillation
)
    : densities_per_m3_(std::move(densities_per_m3)),
      radii_m_(std::move(radii_m)),
      uses_automatic_wavenumber_(true),
      automatic_radial_resolution_m_(automatic_radial_resolution_m),
      automatic_samples_per_oscillation_(automatic_samples_per_oscillation)
{
    require(!densities_per_m3_.empty(), "densities must not be empty.");
    require(densities_per_m3_.size() == radii_m_.size(), "densities and radii must have the same length.");
    require(automatic_radial_resolution_m_ >= 0.0, "radial_resolution must be positive.");
    require(automatic_samples_per_oscillation_ >= 8, "samples_per_oscillation must be at least 8.");

    for (double r : radii_m_) {
        require(std::isfinite(r) && r > 0.0, "radii must be positive and finite.");
    }
    for (double n : densities_per_m3_) {
        require(std::isfinite(n) && n >= 0.0, "densities must be non-negative and finite.");
    }
}

bool PercusYevickSolver::uses_automatic_wavenumber() const {
    return uses_automatic_wavenumber_;
}

std::optional<std::string> PercusYevickSolver::radial_grid_warning(
    const std::vector<double>& distances_m
) const {
    if (wavenumber_per_m_.size() < 2 || distances_m.empty()) {
        return std::nullopt;
    }

    const double maximum_distance = *std::max_element(distances_m.begin(), distances_m.end());
    if (maximum_distance <= 0.0) {
        return std::nullopt;
    }

    double maximum_step = 0.0;
    for (std::size_t index = 1; index < wavenumber_per_m_.size(); ++index) {
        maximum_step = std::max(maximum_step, wavenumber_per_m_[index] - wavenumber_per_m_[index - 1]);
    }

    // The sinc kernel oscillates with period 2*pi/r in reciprocal space.
    // Eight samples per period is a conservative minimum for the trapezoidal
    // inverse transform used by radial_fourier_h_from_H.
    const double samples_per_oscillation = (2.0 * pi_value) / (maximum_step * maximum_distance);
    constexpr double minimum_samples_per_oscillation = 8.0;
    if (samples_per_oscillation >= minimum_samples_per_oscillation) {
        return std::nullopt;
    }

    std::ostringstream message;
    message << "The wavenumber grid is likely too coarse for distances up to "
            << maximum_distance << " m: it provides " << samples_per_oscillation
            << " samples per sinc-kernel oscillation (recommended: at least "
            << minimum_samples_per_oscillation << "). Increase the number of wavenumber "
            << "points or reduce the maximum wavenumber.";
    return message.str();
}

std::vector<double> PercusYevickSolver::compute_epsilons() const {
    std::vector<double> eps(4, 0.0);
    for (int alpha = 0; alpha <= 3; ++alpha) {
        double sum_term = 0.0;
        for (std::size_t i = 0; i < radii_m_.size(); ++i) {
            const double term = densities_per_m3_[i] * std::pow(2.0 * radii_m_[i], static_cast<double>(alpha));
            sum_term += term;
        }
        eps[static_cast<std::size_t>(alpha)] = (pi_value / 6.0) * sum_term;
    }
    return eps;
}

void PercusYevickSolver::compute_parameters(
    const std::vector<double>& epsilons,
    std::vector<double>& R_ij_m,
    std::vector<double>& S_ij_m,
    std::vector<double>& A_i,
    std::vector<double>& B_i_m,
    std::vector<double>& D_ij_m2
) const {
    const std::size_t N = radii_m_.size();
    R_ij_m.assign(N * N, 0.0);
    S_ij_m.assign(N * N, 0.0);

    const double e0 = epsilons[0];
    (void)e0;
    const double e1 = epsilons[1];
    const double e2 = epsilons[2];
    const double e3 = epsilons[3];

    const double denominator = 1.0 - e3;
    require(std::isfinite(denominator) && denominator != 0.0, "Invalid (1 - e3) denominator.");

    A_i.assign(N, 0.0);
    B_i_m.assign(N, 0.0);
    D_ij_m2.assign(N * N, 0.0);

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            const double ri = radii_m_[i];
            const double rj = radii_m_[j];
            R_ij_m[index_2d(i, j, N)] = ri + rj;
            S_ij_m[index_2d(i, j, N)] = rj - ri;
        }
    }

    for (std::size_t i = 0; i < N; ++i) {
        const double r = radii_m_[i];
        A_i[i] = (1.0 - e3 + 6.0 * r * e2) / (denominator * denominator);
        B_i_m[i] = (-6.0 * r * r * e2) / (denominator * denominator);
    }

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            const double Rij = R_ij_m[index_2d(i, j, N)];
            const double Dij = -(A_i[i] * Rij * Rij) / 2.0 - (B_i_m[i] * Rij);
            D_ij_m2[index_2d(i, j, N)] = Dij;
        }
    }
}

std::vector<double> PercusYevickSolver::compute_Cpy(const std::vector<double>& epsilons) const {
    const std::size_t N = radii_m_.size();
    const std::size_t P = wavenumber_per_m_.size();

    const double e0 = epsilons[0];
    const double e1 = epsilons[1];
    const double e2 = epsilons[2];
    const double e3 = epsilons[3];

    const double denom = 1.0 - e3;
    require(std::isfinite(denom) && denom != 0.0, "Invalid (1 - e3) denominator.");

    // Precompute per species and p:
    // X_k(p) = r_k * p
    // R_k = 2 r_k
    // N_k(p) = R_k^2 * sin(X)/X
    // M_k(p) = 3 R_k^3 * (sin/X^3 - cos/X^2)
    std::vector<double> X(N * P, 0.0);
    std::vector<double> Nk(N * P, 0.0);
    std::vector<double> Mk(N * P, 0.0);

    for (std::size_t k = 0; k < N; ++k) {
        const double r = radii_m_[k];
        const double Rk = 2.0 * r;
        for (std::size_t p = 0; p < P; ++p) {
            const double pp = wavenumber_per_m_[p];
            const double x = r * pp;
            X[index_2d(k, p, P)] = x;

            const double s_over_x = safe_sin_over_x(x);
            const double s_over_x3_minus_c_over_x2 = safe_sin_over_x3_minus_cos_over_x2(x);

            Nk[index_2d(k, p, P)] = (Rk * Rk) * s_over_x;
            Mk[index_2d(k, p, P)] = 3.0 * (Rk * Rk * Rk) * s_over_x3_minus_c_over_x2;
        }
    }

    std::vector<double> C(N * N * P, 0.0);

    for (std::size_t i = 0; i < N; ++i) {
        const double ni = densities_per_m3_[i];
        const double Ri = 2.0 * radii_m_[i];

        for (std::size_t j = 0; j < N; ++j) {
            const double nj = densities_per_m3_[j];
            const double Rj = 2.0 * radii_m_[j];

            const double term0 = -(pi_value / 6.0) * std::sqrt(ni * nj) / denom;

            for (std::size_t p = 0; p < P; ++p) {
                const double pp = wavenumber_per_m_[p];

                const double Xi = X[index_2d(i, p, P)];
                const double Xj = X[index_2d(j, p, P)];

                const double Ni_p = Nk[index_2d(i, p, P)];
                const double Nj_p = Nk[index_2d(j, p, P)];

                const double Mi_p = Mk[index_2d(i, p, P)];
                const double Mj_p = Mk[index_2d(j, p, P)];

                const double cos_Xi = std::cos(Xi);
                const double sin_Xi = std::sin(Xi);
                const double cos_Xj = std::cos(Xj);
                const double sin_Xj = std::sin(Xj);

                const double bracket_i =
                    cos_Xi
                    + Xi * sin_Xi
                    + (3.0 * e2 * Ri * cos_Xi) / denom
                    + (3.0 * e1 * Ni_p) / denom
                    + (9.0 * e2 * e2 * Ni_p) / (denom * denom);

                const double bracket_j =
                    cos_Xj
                    + Xj * sin_Xj
                    + (3.0 * e2 * Rj * cos_Xj) / denom
                    + (3.0 * e1 * Nj_p) / denom
                    + (9.0 * e2 * e2 * Nj_p) / (denom * denom);

                const double term2 = Mj_p * bracket_i;
                const double term3 = Mi_p * bracket_j;

                const double term4_factor =
                    (e0 / denom)
                    + (pp * pp * e2) / (4.0 * denom)
                    + (6.0 * e1 * e2) / (denom * denom)
                    + (9.0 * e2 * e2 * e2) / (denom * denom * denom);

                const double term4 = Mi_p * Mj_p * term4_factor;

                const double term5 = 3.0 * Ni_p * Rj * cos_Xj;
                const double term6 = 3.0 * Nj_p * Ri * cos_Xi;

                const double term7 = (9.0 * e2 * Ni_p * Nj_p) / denom;

                const double value = term0 * (term2 + term3 + term4 + term5 + term6 + term7);

                C[index_3d(i, j, p, N, P)] = value;
            }
        }
    }

    return C;
}

std::vector<double> PercusYevickSolver::solve_H_from_C_batch(
    const std::vector<double>& C,
    std::size_t N,
    std::size_t P
) {
    require(C.size() == N * N * P, "C has an invalid size.");

    std::vector<double> H(N * N * P, 0.0);

    // For each p, solve (I - C_p) H_p = C_p
    std::vector<double> A(N * N, 0.0);
    std::vector<double> B(N * N, 0.0);

    for (std::size_t p = 0; p < P; ++p) {
        // Build A and B for this p
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                const double Cij = C[index_3d(i, j, p, N, P)];
                A[index_2d(i, j, N)] = (i == j ? 1.0 : 0.0) - Cij;
                B[index_2d(i, j, N)] = Cij;
            }
        }

        lu_solve_in_place_multiple_rhs(A, B, N);

        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                H[index_3d(i, j, p, N, P)] = B[index_2d(i, j, N)];
            }
        }
    }

    return H;
}

std::vector<double> PercusYevickSolver::radial_fourier_h_from_H(
    const std::vector<double>& H,
    const std::vector<double>& distances_m,
    const std::vector<double>& wavenumber_per_m,
    const std::vector<double>& densities_per_m3,
    std::size_t N
) {
    const std::size_t P = wavenumber_per_m.size();
    const std::size_t R = distances_m.size();

    require(H.size() == N * N * P, "H has an invalid size.");
    require(P >= 2, "wavenumber must have at least two points for trapezoid integration.");

    std::vector<double> h(N * N * R, 0.0);

    // h_ij(r) = (1 / (2 pi^2)) * (1 / sqrt(n_i n_j))
    //           * ∫ H_ij(k) k^2 sin(kr)/(kr) dk
    const double prefactor = 1.0 / (2.0 * pi_value * pi_value);

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            const double ni = densities_per_m3[i];
            const double nj = densities_per_m3[j];

            if (!(ni > 0.0) || !(nj > 0.0)) {
                for (std::size_t r = 0; r < R; ++r) {
                    h[index_3d(i, j, r, N, R)] = 0.0;
                }
                continue;
            }

            const double dens_factor = std::sqrt(ni * nj);


            require(dens_factor > 0.0, "densities must be positive for h(r) conversion.");

            for (std::size_t r = 0; r < R; ++r) {
                const double rr = distances_m[r];

                double integral = 0.0;

                for (std::size_t k = 0; k + 1 < P; ++k) {
                    const double p0 = wavenumber_per_m[k];
                    const double p1 = wavenumber_per_m[k + 1];
                    const double dp = p1 - p0;

                    const double x0 = rr * p0;
                    const double x1 = rr * p1;

                    const double kernel0 = (rr == 0.0) ? 1.0 : safe_sin_over_x(x0);
                    const double kernel1 = (rr == 0.0) ? 1.0 : safe_sin_over_x(x1);

                    const double f0 = H[index_3d(i, j, k, N, P)] * (p0 * p0) * kernel0;
                    const double f1 = H[index_3d(i, j, k + 1, N, P)] * (p1 * p1) * kernel1;

                    integral += 0.5 * dp * (f0 + f1);
                }

                h[index_3d(i, j, r, N, R)] = prefactor * (integral / dens_factor);
            }
        }
    }

    return h;
}

PercusYevickResult PercusYevickSolver::compute(std::vector<double> distances_m) const {
    require(!distances_m.empty(), "distances must not be empty.");

    for (double r : distances_m) {
        require(std::isfinite(r) && r >= 0.0, "distances must be non negative and finite.");
    }

    if (uses_automatic_wavenumber_) {
        const double maximum_distance = *std::max_element(distances_m.begin(), distances_m.end());
        require(maximum_distance > 0.0, "automatic wavenumber selection requires at least one positive distance.");

        const double minimum_radius = *std::min_element(radii_m_.begin(), radii_m_.end());
        const double radial_resolution = automatic_radial_resolution_m_ > 0.0
            ? automatic_radial_resolution_m_
            : minimum_radius / 20.0;
        const double maximum_wavenumber = pi_value / radial_resolution;
        const double wavenumber_step = (2.0 * pi_value) /
            (static_cast<double>(automatic_samples_per_oscillation_) * maximum_distance);
        const std::size_t number_of_points = static_cast<std::size_t>(
            std::ceil(maximum_wavenumber / wavenumber_step)
        ) + 1;

        std::vector<double> automatic_wavenumber(number_of_points, 0.0);
        for (std::size_t index = 0; index < number_of_points; ++index) {
            automatic_wavenumber[index] = maximum_wavenumber *
                static_cast<double>(index) / static_cast<double>(number_of_points - 1);
        }

        return PercusYevickSolver(
            densities_per_m3_, radii_m_, std::move(automatic_wavenumber)
        ).compute(std::move(distances_m));
    }

    const std::size_t N = radii_m_.size();
    const std::size_t P = wavenumber_per_m_.size();
    const std::size_t R = distances_m.size();

    PercusYevickResult result;
    result.number_of_species = N;
    result.number_of_wavenumber_points = P;
    result.number_of_r_points = R;

    result.radii_m = radii_m_;
    result.densities_per_m3 = densities_per_m3_;
    result.wavenumber_per_m = wavenumber_per_m_;
    result.distances_m = std::move(distances_m);

    result.epsilons = compute_epsilons();

    compute_parameters(
        result.epsilons,
        result.R_ij_m,
        result.S_ij_m,
        result.A_i,
        result.B_i_m,
        result.D_ij_m2
    );

    result.Cpy = compute_Cpy(result.epsilons);
    result.H = solve_H_from_C_batch(result.Cpy, N, P);

    result.h = this->radial_fourier_h_from_H(result.H, result.distances_m, result.wavenumber_per_m, result.densities_per_m3, N);

    result.g.assign(N * N * R, 0.0);
    for (std::size_t i = 0; i < N * N * R; ++i) {
        result.g[i] = 1.0 + result.h[i];
    }

    return result;
}
