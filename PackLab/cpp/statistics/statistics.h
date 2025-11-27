#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>

struct Statistics {
    std::size_t attempted_insertions = 0;
    std::size_t accepted_insertions = 0;
    std::size_t rejected_insertions = 0;
    std::size_t consecutive_rejections = 0;

    std::size_t sphere_count = 0;

    double packing_fraction_geometry = 0.0;
    double packing_fraction_simulator = 0.0;

    double radius_min = 0.0;
    double radius_max = 0.0;
    double radius_mean = 0.0;
    double radius_median = 0.0;
    double radius_std = 0.0;

    double total_runtime_seconds = 0.0;

    void print() const;
};