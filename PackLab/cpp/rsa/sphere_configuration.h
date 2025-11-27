class SphereConfiguration {
public:
    const std::vector<Vector3d>& center_positions() const { return center_positions_values_; }
    const std::vector<double>& radii() const { return radii_values_; }

    double total_sphere_volume() const {
        double volume_sum = 0.0;
        for (double radius : radii_values_) {
            volume_sum += (4.0 / 3.0) * PI * radius * radius * radius;
        }
        return volume_sum;
    }

private:
    friend class Simulator;

    std::vector<Vector3d> center_positions_values_;
    std::vector<double> radii_values_;
};





