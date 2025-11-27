#include "domain.h"

Vector3d Domain::wrap_position_if_periodic(const Vector3d& position) const {
    if (!use_periodic_boundaries_value_)
        return position;


    auto wrap_01 = [](double value, double length) {
        double wrapped = std::fmod(value, length);
        if (wrapped < 0.0) wrapped += length;
        return wrapped;
    };

    return Vector3d{
        wrap_01(position.x, length_x_value_),
        wrap_01(position.y, length_y_value_),
        wrap_01(position.z, length_z_value_)
    };
}

Vector3d Domain::sample_uniform_position(std::mt19937_64& random_generator, double margin) const {
    const double effective_margin = std::max(0.0, margin);

    const double minimum_x = use_periodic_boundaries_value_ ? 0.0 : effective_margin;
    const double minimum_y = use_periodic_boundaries_value_ ? 0.0 : effective_margin;
    const double minimum_z = use_periodic_boundaries_value_ ? 0.0 : effective_margin;

    const double maximum_x = use_periodic_boundaries_value_ ? length_x_value_ : (length_x_value_ - effective_margin);
    const double maximum_y = use_periodic_boundaries_value_ ? length_y_value_ : (length_y_value_ - effective_margin);
    const double maximum_z = use_periodic_boundaries_value_ ? length_z_value_ : (length_z_value_ - effective_margin);

    if (!use_periodic_boundaries_value_) {
        if (maximum_x <= minimum_x || maximum_y <= minimum_y || maximum_z <= minimum_z) {
            throw std::invalid_argument("Margin is too large for the domain size.");
        }
    }

    std::uniform_real_distribution<double> uniform_x(minimum_x, maximum_x);
    std::uniform_real_distribution<double> uniform_y(minimum_y, maximum_y);
    std::uniform_real_distribution<double> uniform_z(minimum_z, maximum_z);

    return wrap_position_if_periodic(Vector3d{uniform_x(random_generator), uniform_y(random_generator), uniform_z(random_generator)});
}

double Domain::minimum_image_displacement(double delta, double box_length) const {
    if (!use_periodic_boundaries_value_)
        return delta;

    const double half_length = 0.5 * box_length;

    if (delta > half_length)
        return delta - box_length;
    if (delta < -half_length)
        return delta + box_length;
    return delta;
}

