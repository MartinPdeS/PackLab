// domain_box.hpp
#pragma once

#include <vector> // for std::vector
#include <stdexcept> // for std::invalid_argument
#include <random> // for std::mt19937_64

#include "../utils/utils.h"


class Domain {
public:
    Domain() = default;

    Domain(double _length_x, double _length_y, double _length_z, bool _use_periodic_boundaries)
    :   length_x(_length_x),
        length_y(_length_y),
        length_z(_length_z),
        use_periodic_boundaries(_use_periodic_boundaries)
    {
        if (length_x <= 0.0 || length_y <= 0.0 || length_z <= 0.0)
            throw std::invalid_argument("Box lengths must be positive.");

        this->volume = this->get_volume();
    }

    double get_volume() const {return length_x * length_y * length_z;}

    Vector3d wrap_position_if_periodic(const Vector3d& position) const;

    Vector3d sample_uniform_position(std::mt19937_64& random_generator, double margin) const;

    double minimum_image_displacement(double delta, double box_length) const;

    void scale( double scale_factor ) {
        if (scale_factor <= 0.0)
            throw std::invalid_argument("scale_factor must be positive.");

        length_x *= scale_factor;
        length_y *= scale_factor;
        length_z *= scale_factor;
        volume = this->get_volume();
    }

    double length_x = 1.0;
    double length_y = 1.0;
    double length_z = 1.0;
    double volume = 1.0;
    bool use_periodic_boundaries = true;
};
