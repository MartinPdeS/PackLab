#pragma once

#include <cstdint>    // for std::int64_t
#include <cstddef>    // for std::size_t, std::int64_t
#include <functional> // for std::hash
#include <cmath>      // for std::sqrt, std::exp

static const double PI = 3.141592653589793238462643383279502884;

struct Vector3d {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vector3d operator-(const Vector3d& other) const {
        return Vector3d{ x - other.x, y - other.y, z - other.z };
    }

    Vector3d operator+(const Vector3d& other) const {
        return Vector3d{ x + other.x, y + other.y, z + other.z };
    }

    double norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }
};


struct CellIndex {
    std::int64_t i = 0;
    std::int64_t j = 0;
    std::int64_t k = 0;
};


struct CellIndexHasher {
    std::size_t operator()(const CellIndex& cell_index) const noexcept {
        // Basic mixing for three int64 values
        std::size_t h1 = std::hash<std::int64_t>{}(cell_index.i);
        std::size_t h2 = std::hash<std::int64_t>{}(cell_index.j);
        std::size_t h3 = std::hash<std::int64_t>{}(cell_index.k);
        std::size_t mixed = h1;
        mixed ^= h2 + 0x9e3779b97f4a7c15ULL + (mixed << 6) + (mixed >> 2);
        mixed ^= h3 + 0x9e3779b97f4a7c15ULL + (mixed << 6) + (mixed >> 2);
        return mixed;
    }
};

struct CellIndexEqual {
    bool operator()(const CellIndex& a, const CellIndex& b) const noexcept {
        return a.i == b.i && a.j == b.j && a.k == b.k;
    }
};


struct Options {
    std::uint64_t random_seed = 0;

    std::size_t maximum_attempts = 2'000'000;
    std::size_t maximum_spheres = 0; // 0 means unlimited

    std::size_t maximum_consecutive_rejections = 200'000;

    double target_packing_fraction = 0.0; // 0 means ignore
    double minimum_center_separation_addition = 0.0; // extra added to (r_i + r_j)
    double containment_padding = 0.0; // extra margin inside walls, ignored for periodic

    double spatial_grid_cell_size = 0.0; // 0 means auto based on radii seen so far

    bool enforce_radii_distribution = true;
};
