// domain_box.hpp
#pragma once

#include <memory>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "../utils/utils.h"


class Domain {
public:
    Domain() = default;

    Domain(double length_x, double length_y, double length_z, bool use_periodic_boundaries)
    :   length_x_value_(length_x),
        length_y_value_(length_y),
        length_z_value_(length_z),
        use_periodic_boundaries_value_(use_periodic_boundaries)
    {
        if (length_x_value_ <= 0.0 || length_y_value_ <= 0.0 || length_z_value_ <= 0.0)
            throw std::invalid_argument("Box lengths must be positive.");
    }

    double length_x() const { return length_x_value_; }
    double length_y() const { return length_y_value_; }
    double length_z() const { return length_z_value_; }
    bool use_periodic_boundaries() const { return use_periodic_boundaries_value_; }
    double volume() const {return length_x_value_ * length_y_value_ * length_z_value_;}

    Vector3d wrap_position_if_periodic(const Vector3d& position) const;

    Vector3d sample_uniform_position(std::mt19937_64& random_generator, double margin) const;

    double minimum_image_displacement(double delta, double box_length) const;

private:
    double length_x_value_ = 1.0;
    double length_y_value_ = 1.0;
    double length_z_value_ = 1.0;
    bool use_periodic_boundaries_value_ = true;
};
