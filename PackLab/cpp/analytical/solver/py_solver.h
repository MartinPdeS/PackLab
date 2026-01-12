#pragma once

#include <cstddef>
#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>


struct PercusYevickResult {
    std::vector<double> epsilons;               // size 4

    std::vector<double> radii_m;                // size N
    std::vector<double> densities_per_m3;       // size N

    std::vector<double> p_per_m;                // size P
    std::vector<double> distances_m;            // size R

    std::vector<double> R_ij_m;                 // size N*N
    std::vector<double> S_ij_m;                 // size N*N
    std::vector<double> A_i;                    // size N
    std::vector<double> B_i_m;                  // size N
    std::vector<double> D_ij_m2;                // size N*N

    std::vector<double> Cpy;                    // size N*N*P
    std::vector<double> H;                      // size N*N*P

    std::vector<double> h;                      // size N*N*R
    std::vector<double> g;                      // size N*N*R

    std::size_t number_of_species = 0;          // N
    std::size_t number_of_p_points = 0;         // P
    std::size_t number_of_r_points = 0;         // R
};

class PercusYevickSolver {
public:
    PercusYevickSolver(
        std::vector<double> densities_per_m3,
        std::vector<double> radii_m,
        std::vector<double> p_per_m
    );

    PercusYevickResult compute(std::vector<double> distances_m) const;

private:
    std::vector<double> compute_epsilons() const;
    void compute_parameters(
        const std::vector<double>& epsilons,
        std::vector<double>& R_ij_m,
        std::vector<double>& S_ij_m,
        std::vector<double>& A_i,
        std::vector<double>& B_i_m,
        std::vector<double>& D_ij_m2
    ) const;

    std::vector<double> compute_Cpy(const std::vector<double>& epsilons) const;
    static std::vector<double> solve_H_from_C_batch(
        const std::vector<double>& C,
        std::size_t N,
        std::size_t P
    );

    static std::vector<double> radial_fourier_h_from_H(
        const std::vector<double>& H,
        const std::vector<double>& distances_m,
        const std::vector<double>& p_per_m,
        const std::vector<double>& densities_per_m3,
        std::size_t N
    );

private:
    std::vector<double> densities_per_m3_;
    std::vector<double> radii_m_;
    std::vector<double> p_per_m_;
};
