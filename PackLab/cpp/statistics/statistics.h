#pragma once

#include <cstddef>
#include <iostream>
#include <iomanip>

struct Statistics {
    // Counters
    std::size_t attempted_insertions = 0;
    std::size_t accepted_insertions = 0;
    std::size_t rejected_insertions = 0;
    std::size_t consecutive_rejections = 0;

    std::size_t sphere_count = 0;

    // Radius information
    double radius_min = 0.0;
    double radius_max = 0.0;
    double radius_mean = 0.0;
    double radius_median = 0.0;
    double radius_std = 0.0;

    // Packing fractions
    double packing_fraction_geometry = 0.0;
    double packing_fraction_simulator = 0.0;

    // Runtime
    double total_runtime_seconds = 0.0;

private:
    // Internal timestamps stored as double (seconds, high resolution)
    double benchmark_start_seconds = 0.0;
    double benchmark_end_seconds = 0.0;

public:
    void start_benchmark();
    void end_benchmark();
    void print() const;
};
