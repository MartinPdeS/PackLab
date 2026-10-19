#pragma once


class SphereConfiguration {
public:
    std::vector<Vector3d> center_positions;
    std::vector<double> radii_values;
    std::vector<int> class_index_values;

    const std::vector<double>& radii() const { return radii_values; }

    double total_sphere_volume() const {
        double volume_sum = 0.0;
        for (double radius : radii_values) {
            volume_sum += (4.0 / 3.0) * PI * radius * radius * radius;
        }
        return volume_sum;
    }

    friend class Simulator;

};

